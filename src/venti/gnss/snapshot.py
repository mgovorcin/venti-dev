# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Frozen UNR grid snapshots (PRD R-G4, R-O4; plan T48).

Operations calibrate against a *frozen* copy of the UNR gridded GNSS
product, not the live server: the same inputs, snapshot and config must give
the same product (R-O4), and a new UNR release rolls forward deliberately,
through the e2e gate (plan T46), roughly every six months.

A snapshot is a directory::

    unr_grid_<version>_<YYYYMMDD>/
        snapshot.json            what it is: id, UNR version, frame, grid types,
                                 node counts, data span, hashes, tool versions
        MANIFEST.sha256          one line per file, ``sha256sum -c`` compatible
        grid_latlon_lookup.txt   the UNR lookup, byte for byte
        nodes/<id:06d>_<frame>_<grid_type>.tenv8

``nodes/`` is exactly the layout `venti.gnss.sampling.GnssGridConfig`
reads, so `grid_config_from_snapshot` is all a workflow needs, and the
snapshot id and lookup hash flow into the GNSS `Provenance` of every
product. Retrieval goes through geepers (plan D10: single source for UNR
access).
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import shutil
import tempfile
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

import numpy as np

from .sampling import GnssGridConfig

logger = logging.getLogger(__name__)

__all__ = [
    "LookupFetcher",
    "NodeDownloader",
    "SnapshotInfo",
    "grid_config_from_snapshot",
    "load_snapshot",
    "snapshot_unr_grid",
    "verify_snapshot",
]

LOOKUP_NAME = "grid_latlon_lookup.txt"
MANIFEST_NAME = "MANIFEST.sha256"
INFO_NAME = "snapshot.json"
NODES_DIR = "nodes"

# ``fetch(url) -> text`` and ``download(ids, out_dir, grid_type) -> files``
LookupFetcher = Callable[[str], str]
NodeDownloader = Callable[[Sequence[int], Path, str], list[Path]]


@dataclass
class SnapshotInfo:
    """Contents of ``snapshot.json``."""

    snapshot_id: str
    unr_version: str
    reference_frame: str
    grid_types: list[str]
    created: str
    lookup_sha256: str
    n_lookup_nodes: int
    n_selected_nodes: int
    n_downloaded: dict[str, int]
    n_failed: dict[str, int]
    failed_ids: dict[str, list[int]]
    data_span: dict[str, list[float] | None]
    bounds_snwe: list[float] | None
    source: dict[str, str]
    tool_versions: dict[str, str]
    layout: str = "nodes/<id:06d>_<frame>_<grid_type>.tenv8"
    notes: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class _Selection:
    ids: np.ndarray
    n_lookup: int
    bounds: list[float] | None = None
    extra: dict[str, Any] = field(default_factory=dict)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _pkg_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _geepers_urls(version: str) -> dict[str, str]:
    from geepers.gps_sources import unr_grid as g

    return {
        "lookup_url": g.LOOKUP_FILE_URL.format(version=version),
        "data_url_template": g.GRID_DATA_BASE_URL,
        "gridded_dirs": json.dumps(g.GRIDDED_TYPE_DIRS),
    }


def _default_fetch_lookup(url: str) -> str:
    import requests

    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    return resp.text


def _default_downloader(
    version: str, reference_frame: str, max_workers: int
) -> NodeDownloader:
    """Download through geepers; a failing node is retried alone and recorded."""

    def download(ids: Sequence[int], out_dir: Path, grid_type: str) -> list[Path]:
        from geepers.gps_sources.unr_grid import UnrGridSource

        src = UnrGridSource(version=version, gridded_type=grid_type, cache_dir=out_dir)  # type: ignore[arg-type]
        id_strs = [f"{int(i):06d}" for i in ids]
        kwargs: dict[str, Any] = {
            "plate": reference_frame,
            "output_dir": out_dir,
            "version": version,
            "gridded_type": grid_type,
        }
        try:
            return list(
                src.download_data_files(id_strs, max_workers=max_workers, **kwargs)
            )
        except Exception as exc:
            logger.warning("bulk download failed (%s); retrying node by node", exc)
        files: list[Path] = []
        for s in id_strs:
            try:
                files.extend(src.download_data_files([s], max_workers=1, **kwargs))
            except Exception as exc:
                logger.warning("node %s: %s", s, exc)
        return files

    return download


def _parse_lookup(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = [ln.split() for ln in text.splitlines() if ln.strip()]
    if not rows:
        msg = "the grid lookup is empty"
        raise ValueError(msg)
    arr = np.array([[float(r[0]), float(r[1]), float(r[2])] for r in rows])
    ids = arr[:, 0].astype(np.int64)
    lon = ((arr[:, 1] + 180.0) % 360.0) - 180.0  # v0.3 publishes 0..360
    return ids, lon, arr[:, 2]


def _select(text: str, bounds_snwe: Sequence[float] | None) -> _Selection:
    ids, lon, lat = _parse_lookup(text)
    if bounds_snwe is None:
        return _Selection(ids=ids, n_lookup=len(ids))
    s, n, w, e = (float(v) for v in bounds_snwe)
    if not (s < n and w < e):
        msg = f"bounds must be (S, N, W, E) with S < N and W < E, got {bounds_snwe}"
        raise ValueError(msg)
    keep = (lat >= s) & (lat <= n) & (lon >= w) & (lon <= e)
    return _Selection(ids=ids[keep], n_lookup=len(ids), bounds=[s, n, w, e])


def _data_span(files: Sequence[Path]) -> list[float] | None:
    lo, hi = np.inf, -np.inf
    for f in files:
        with f.open() as fh:
            first = last = ""
            for line in fh:
                if line.strip():
                    if not first:
                        first = line
                    last = line
        if not first:
            continue
        lo = min(lo, float(first.split()[0]))
        hi = max(hi, float(last.split()[0]))
    return None if not np.isfinite(lo) else [lo, hi]


def _write_manifest(snapshot_dir: Path) -> int:
    lines = []
    for path in sorted(p for p in snapshot_dir.rglob("*") if p.is_file()):
        if path.name == MANIFEST_NAME:
            continue
        lines.append(
            f"{_sha256_file(path)}  {path.relative_to(snapshot_dir).as_posix()}"
        )
    (snapshot_dir / MANIFEST_NAME).write_text("\n".join(lines) + "\n")
    return len(lines)


def snapshot_unr_grid(
    output_root: Path | str,
    *,
    version: str = "0.3",
    reference_frame: str = "IGS20",
    grid_types: Sequence[str] = ("constant",),
    bounds_snwe: Sequence[float] | None = None,
    snapshot_date: date | None = None,
    max_workers: int = 8,
    notes: str = "",
    fetch_lookup: LookupFetcher | None = None,
    download: NodeDownloader | None = None,
) -> Path:
    """Write a frozen snapshot of the UNR grid and return its directory.

    Parameters
    ----------
    output_root : Path
        Parent directory; the snapshot is ``unr_grid_<version>_<YYYYMMDD>``.
    version, reference_frame : str
        UNR product version (``'0.3'``) and frame (``'IGS20'``; R-G1).
    grid_types : sequence of str
        ``'constant'`` (operations) and/or ``'variable'`` (reprocessing).
    bounds_snwe : (S, N, W, E), optional
        Only nodes inside these degrees; default the whole grid.
    snapshot_date : date, optional
        Date in the snapshot id; default today (UTC).
    max_workers : int
        Parallel downloads.
    notes : str
        Free text recorded in ``snapshot.json``.
    fetch_lookup, download : callables, optional
        Replace the network access (tests, mirrors).

    Raises
    ------
    FileExistsError
        If the snapshot directory already exists (snapshots are immutable).

    """
    root = Path(output_root)
    day = snapshot_date or datetime.now(UTC).date()
    snapshot_id = f"unr_grid_{version}_{reference_frame}_{day:%Y%m%d}"
    snapshot_dir = root / f"unr_grid_{version}_{day:%Y%m%d}"
    if snapshot_dir.exists():
        msg = f"{snapshot_dir} exists; a snapshot is never overwritten"
        raise FileExistsError(msg)
    for gt in grid_types:
        if gt not in ("constant", "variable"):
            msg = f"unknown grid type {gt!r}"
            raise ValueError(msg)
    if bounds_snwe is not None:
        _select("0 0 0", bounds_snwe)  # validates the bounds before anything is written
    urls = _geepers_urls(version)
    fetch = fetch_lookup or _default_fetch_lookup
    get_nodes = download or _default_downloader(version, reference_frame, max_workers)

    snapshot_dir.mkdir(parents=True)
    try:
        _build(
            snapshot_dir,
            snapshot_id,
            version,
            reference_frame,
            grid_types,
            bounds_snwe,
            urls,
            fetch,
            get_nodes,
            notes,
        )
    except BaseException:
        # a half-written directory must never pass for a snapshot
        shutil.rmtree(snapshot_dir, ignore_errors=True)
        raise
    return snapshot_dir


def _build(
    snapshot_dir: Path,
    snapshot_id: str,
    version: str,
    reference_frame: str,
    grid_types: Sequence[str],
    bounds_snwe: Sequence[float] | None,
    urls: dict[str, str],
    fetch: LookupFetcher,
    get_nodes: NodeDownloader,
    notes: str,
) -> None:
    nodes_dir = snapshot_dir / NODES_DIR
    nodes_dir.mkdir()
    lookup_text = fetch(urls["lookup_url"])
    lookup_path = snapshot_dir / LOOKUP_NAME
    lookup_path.write_text(lookup_text, encoding="utf-8")
    sel = _select(lookup_text, bounds_snwe)
    logger.info(
        "snapshot %s: %d of %d lookup nodes selected",
        snapshot_id,
        len(sel.ids),
        sel.n_lookup,
    )

    n_downloaded: dict[str, int] = {}
    n_failed: dict[str, int] = {}
    failed_ids: dict[str, list[int]] = {}
    span: dict[str, list[float] | None] = {}
    for gt in grid_types:
        with tempfile.TemporaryDirectory(dir=snapshot_dir, prefix=f".dl_{gt}_") as tmp:
            files = get_nodes([int(i) for i in sel.ids], Path(tmp), gt)
            got: set[int] = set()
            kept: list[Path] = []
            for f in files:
                node = int(Path(f).name.split("_")[0])
                dest = nodes_dir / f"{node:06d}_{reference_frame}_{gt}.tenv8"
                shutil.move(str(f), dest)
                got.add(node)
                kept.append(dest)
        missing = sorted(int(i) for i in sel.ids if int(i) not in got)
        n_downloaded[gt] = len(kept)
        n_failed[gt] = len(missing)
        failed_ids[gt] = missing[:1000]
        span[gt] = _data_span(kept)
        logger.info(
            "%s: %d nodes written, %d failed, data span %s",
            gt,
            len(kept),
            len(missing),
            span[gt],
        )

    info = SnapshotInfo(
        snapshot_id=snapshot_id,
        unr_version=version,
        reference_frame=reference_frame,
        grid_types=list(grid_types),
        created=datetime.now(UTC).isoformat(timespec="seconds"),
        lookup_sha256=_sha256_file(lookup_path),
        n_lookup_nodes=sel.n_lookup,
        n_selected_nodes=len(sel.ids),
        n_downloaded=n_downloaded,
        n_failed=n_failed,
        failed_ids=failed_ids,
        data_span=span,
        bounds_snwe=sel.bounds,
        source=urls,
        tool_versions={
            "venti": _pkg_version("venti"),
            "geepers": _pkg_version("geepers"),
        },
        notes=notes,
    )
    (snapshot_dir / INFO_NAME).write_text(json.dumps(info.as_dict(), indent=1) + "\n")
    n_files = _write_manifest(snapshot_dir)
    logger.info("%s: %d files in %s", snapshot_id, n_files, MANIFEST_NAME)


def load_snapshot(snapshot_dir: Path | str) -> SnapshotInfo:
    """Read ``snapshot.json``."""
    with (Path(snapshot_dir) / INFO_NAME).open() as fh:
        return SnapshotInfo(**json.load(fh))


def verify_snapshot(snapshot_dir: Path | str) -> list[str]:
    """Check every file against ``MANIFEST.sha256``; return the problems, [] if OK."""
    root = Path(snapshot_dir)
    manifest = root / MANIFEST_NAME
    if not manifest.exists():
        return [f"missing {MANIFEST_NAME}"]
    problems: list[str] = []
    listed: set[str] = set()
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        digest, rel = line.split(None, 1)
        rel = rel.strip()
        listed.add(rel)
        path = root / rel
        if not path.exists():
            problems.append(f"missing {rel}")
        elif _sha256_file(path) != digest:
            problems.append(f"changed {rel}")
    for path in root.rglob("*"):
        if path.is_file() and path.name != MANIFEST_NAME:
            rel = path.relative_to(root).as_posix()
            if rel not in listed:
                problems.append(f"unlisted {rel}")
    return problems


def grid_config_from_snapshot(
    snapshot_dir: Path | str,
    utm_epsg: int,
    grid_type: str = "constant",
    *,
    verify: bool = False,
    **options: Any,
) -> GnssGridConfig:
    """Build the `GnssGridConfig` that samples this snapshot.

    `options` are the sampling choices (`buffer_meters`,
    `exclude_defo_nodes`, ...). The lookup, node directory, UNR version,
    frame and `snapshot_id` come from the snapshot. With ``verify=True`` the
    manifest is checked first and any problem raises.
    """
    root = Path(snapshot_dir)
    info = load_snapshot(root)
    if grid_type not in info.grid_types:
        msg = f"snapshot {info.snapshot_id} holds {info.grid_types}, not {grid_type!r}"
        raise ValueError(msg)
    if verify:
        problems = verify_snapshot(root)
        if problems:
            msg = f"snapshot {info.snapshot_id} fails verification: {problems[:5]}"
            raise ValueError(msg)
    return GnssGridConfig(
        grid_lookup=root / LOOKUP_NAME,
        station_dir=root / NODES_DIR,
        utm_epsg=utm_epsg,
        reference_frame=info.reference_frame,
        grid_type=grid_type,  # type: ignore[arg-type]
        version=info.unr_version,
        snapshot_id=info.snapshot_id,
        **options,
    )
