# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""The single GNSS sampling path for DISP-CAL and VLM (plan T28, PRD §4.3).

`sample_gnss_enu` turns the frozen UNR grid snapshot into east/north/up
velocity (and sigma) fields on the DISP grid:

1. select grid nodes inside the frame bounds **plus a buffer** (R-G2), so the
   field near the frame edges is constrained by real nodes instead of being
   an extrapolation that moves whenever an interior node changes;
2. read each node's rate from its constant-grid tenv8 file (R-G1; the
   variable grid is allowed only in reprocessing mode);
3. optionally drop the nodes inside curated defo/event areas and
   re-estimate them from the surrounding nodes with GPS Imaging (R-G5,
   `geepers.gps_imaging.reinterpolate_nodes`);
4. interpolate E, N, U and their sigmas onto the target grid;
5. record provenance (grid version/type/frame, snapshot id, lookup hash,
   node counts, a hash of the sampled field) so DISP-CAL can write it to the
   product metadata and VLM can check it (R-G4, §7.3).

The LOS projection happens afterwards, per pixel, with the LOS unit-vector
rasters (`project_field_to_los`). Because E/N/U are gridded first, nodes
outside the swath never need a LOS look, which makes the trade-study LOS
extrapolation unnecessary.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from ..spatial.interpolation import (
    GridLike,
    RbfFunction,
    grid_coordinates,
    interpolate_rbf,
)
from .unr import calculate_station_velocity, find_stations_in_bounds, read_station_rate

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "GnssField",
    "GnssGridConfig",
    "Provenance",
    "project_field_to_los",
    "sample_gnss_enu",
]

VALUE_COLUMNS: tuple[tuple[str, str], ...] = (
    ("ve", "sigma_ve"),
    ("vn", "sigma_vn"),
    ("vu", "sigma_vu"),
)


@dataclass(frozen=True)
class GnssGridConfig:
    """Where the frozen UNR grid snapshot is and how to sample it.

    Attributes
    ----------
    grid_lookup : Path
        UNR ``grid_latlon_lookup.txt`` of the snapshot.
    station_dir : Path
        Directory with the node files ``<id:06d>_<frame>_<grid_type>.tenv8``
        (the layout `venti.gnss.unr.download_station` writes).
    utm_epsg : int
        EPSG of the target grid (the DISP frame's UTM zone).
    reference_frame, grid_type, version
        UNR product coordinates (R-G1). ``'variable'`` needs `reprocessing`.
    buffer_meters : float
        Nodes this far outside the grid bounds are used too (R-G2).
    exclude_defo_nodes, reinterpolate_excluded : bool
        GNSS-side remove-restore (R-G5).
    reprocessing : bool
        Allows ``grid_type='variable'``.
    snapshot_id : str or None
        Identifier of the frozen snapshot (R-G4), recorded in provenance.
    start_year : float
        Earliest epoch used when a rate is fitted (variable grid only).
    rbf_function : str
        Basis of the node-to-grid interpolation.

    """

    grid_lookup: Path
    station_dir: Path
    utm_epsg: int
    reference_frame: str = "IGS20"
    grid_type: Literal["constant", "variable"] = "constant"
    version: str = "0.3"
    buffer_meters: float = 0.0
    exclude_defo_nodes: bool = False
    reinterpolate_excluded: bool = False
    reprocessing: bool = False
    snapshot_id: str | None = None
    start_year: float = 2014.0
    rbf_function: RbfFunction = "cubic"

    def station_file(self, node_id: int) -> Path:
        """Path of a node's tenv8 file in `station_dir`."""
        return (
            Path(self.station_dir)
            / f"{node_id:06d}_{self.reference_frame}_{self.grid_type}.tenv8"
        )

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form (paths as strings)."""
        d = asdict(self)
        d["grid_lookup"] = str(self.grid_lookup)
        d["station_dir"] = str(self.station_dir)
        return d


@dataclass(frozen=True)
class Provenance:
    """What went into a sampled GNSS field; written to product metadata."""

    grid_version: str
    grid_type: str
    reference_frame: str
    snapshot_id: str | None
    lookup_sha256: str
    buffer_meters: float
    n_nodes: int
    n_missing: int
    n_excluded: int
    n_reinterpolated: int
    defo_db_version: str | None
    field_sha256: str
    config_sha256: str

    @property
    def digest(self) -> str:
        """One hash for the whole realisation: config and field together."""
        return hashlib.sha256(
            (self.config_sha256 + self.field_sha256).encode()
        ).hexdigest()

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form including `digest`."""
        return {**asdict(self), "digest": self.digest}


@dataclass
class GnssField:
    """E/N/U velocity and sigma fields (mm/yr) on the target grid, plus the nodes."""

    ve: np.ndarray
    vn: np.ndarray
    vu: np.ndarray
    sigma_ve: np.ndarray
    sigma_vn: np.ndarray
    sigma_vu: np.ndarray
    nodes: pd.DataFrame
    provenance: Provenance
    grid_shape: tuple[int, int] = field(init=False)

    def __post_init__(self) -> None:
        self.grid_shape = self.ve.shape
        for name in ("vn", "vu", "sigma_ve", "sigma_vn", "sigma_vu"):
            if getattr(self, name).shape != self.grid_shape:
                msg = f"{name} shape {getattr(self, name).shape} != {self.grid_shape}"
                raise ValueError(msg)


# --------------------------------------------------------------------------


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_arrays(*arrays: np.ndarray) -> str:
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a, dtype=np.float32).tobytes())
    return h.hexdigest()


def _bounds_snwe(grid: GridLike) -> tuple[float, float, float, float]:
    gx, gy = grid_coordinates(grid)
    gx = np.asarray(gx, float)
    gy = np.asarray(gy, float)
    return (float(gy.min()), float(gy.max()), float(gx.min()), float(gx.max()))


def _read_rates(cfg: GnssGridConfig, node_ids: np.ndarray) -> tuple[pd.DataFrame, int]:
    """Rates per node from the snapshot files; nodes without a file are skipped."""
    import pandas as pd

    rows: list[dict[str, float]] = []
    missing = 0
    for node_id in node_ids:
        path = cfg.station_file(int(node_id))
        if not path.exists():
            missing += 1
            continue
        if cfg.grid_type == "constant":
            ve, vn, vu, se, sn, su = read_station_rate(path)
        else:
            ve, vn, vu, se, sn, su = calculate_station_velocity(path, cfg.start_year)
        rows.append(
            {
                "id": int(node_id),
                "ve": ve,
                "vn": vn,
                "vu": vu,
                "sigma_ve": se,
                "sigma_vn": sn,
                "sigma_vu": su,
            }
        )
    if missing:
        logger.warning(
            "%d of %d grid nodes have no %s file in %s and are skipped",
            missing,
            len(node_ids),
            cfg.grid_type,
            cfg.station_dir,
        )
    return pd.DataFrame(rows), missing


def sample_gnss_enu(
    cfg: GnssGridConfig,
    grid: GridLike,
    exclude: Any | None = None,
    defo_db_version: str | None = None,
) -> GnssField:
    """Sample the UNR grid snapshot as E/N/U velocity fields on `grid`.

    Parameters
    ----------
    cfg : GnssGridConfig
        Snapshot location and sampling options.
    grid : GridLike
        Target grid (GeoTIFF/NetCDF path, dataset, or ``(x, y)`` coordinate
        arrays) in the UTM CRS `cfg.utm_epsg`.
    exclude : shapely geometry, iterable of geometries or GeoSeries, optional
        Defo/event areas in EPSG:4326 (R-G5). Used only when
        `cfg.exclude_defo_nodes` is True.
    defo_db_version : str, optional
        Version of the defo/event database `exclude` came from (provenance).

    Returns
    -------
    GnssField
        Fields in mm/yr with the nodes used and the provenance record.

    Raises
    ------
    ValueError
        ``grid_type='variable'`` outside reprocessing mode (R-G1), or no node
        with data inside the (buffered) bounds.

    """
    if cfg.grid_type == "variable" and not cfg.reprocessing:
        msg = (
            "grid_type='variable' carries coseismic, postseismic, hydrological and "
            "instrumental signal and is allowed only with reprocessing=True (PRD R-G1)"
        )
        raise ValueError(msg)

    bounds = _bounds_snwe(grid)
    nodes_gdf = find_stations_in_bounds(
        Path(cfg.grid_lookup), bounds, cfg.utm_epsg, padding=cfg.buffer_meters
    )
    if nodes_gdf.empty:
        msg = (
            f"No UNR grid nodes within {bounds} (+{cfg.buffer_meters} m) in"
            f" {cfg.grid_lookup}"
        )
        raise ValueError(msg)

    rates, n_missing = _read_rates(cfg, nodes_gdf.index.to_numpy())
    if rates.empty:
        msg = (
            f"None of the {len(nodes_gdf)} nodes in bounds has a file in"
            f" {cfg.station_dir}"
        )
        raise ValueError(msg)
    # the join returns a plain DataFrame, so take the UTM coordinates first
    coords = nodes_gdf[["lon", "lat"]].copy()
    coords["x"] = nodes_gdf.geometry.x.to_numpy()
    coords["y"] = nodes_gdf.geometry.y.to_numpy()
    nodes = rates.set_index("id").join(coords, how="inner").reset_index()

    n_excluded = n_reinterpolated = 0
    nodes["excluded"] = False
    nodes["reinterpolated"] = False
    if cfg.exclude_defo_nodes and exclude is not None:
        from geepers._optional import require

        shapely = require("shapely")
        if hasattr(exclude, "union_all"):
            geom = exclude.union_all()
        elif hasattr(exclude, "geometry"):
            geom = exclude.geometry.union_all()
        elif isinstance(exclude, shapely.Geometry):
            geom = exclude
        else:
            geom = shapely.union_all(list(exclude))
        inside = np.asarray(
            shapely.contains_xy(geom, nodes["lon"].to_numpy(), nodes["lat"].to_numpy()),
            dtype=bool,
        )
        n_excluded = int(inside.sum())
        nodes["excluded"] = inside
        if n_excluded and cfg.reinterpolate_excluded:
            from geepers.gps_imaging import reinterpolate_nodes

            nodes = reinterpolate_nodes(
                nodes, geom, columns=VALUE_COLUMNS, lon="lon", lat="lat"
            )
            n_reinterpolated = int(nodes["reinterpolated"].sum())
        elif n_excluded:
            nodes = nodes[~inside].reset_index(drop=True)
        if nodes.empty:
            msg = "every grid node in bounds lies inside the exclusion areas"
            raise ValueError(msg)
        logger.info(
            "GNSS-side remove-restore: %d nodes inside exclusion areas, %d"
            " re-interpolated, %d dropped",
            n_excluded,
            n_reinterpolated,
            n_excluded - n_reinterpolated,
        )

    x = nodes["x"].to_numpy(float)
    y = nodes["y"].to_numpy(float)
    fields = {
        col: interpolate_rbf(
            grid, x, y, nodes[col].to_numpy(float), function=cfg.rbf_function
        )
        for pair in VALUE_COLUMNS
        for col in pair
    }
    for pair in VALUE_COLUMNS:
        # a gridded sigma is still a sigma: never negative
        fields[pair[1]] = np.maximum(fields[pair[1]], 0.0)

    # portable: the file paths differ between machines, the lookup's content
    # is covered by lookup_sha256 (VLM compares this hash with DISP-CAL's, T55)
    portable = {
        k: v
        for k, v in cfg.as_dict().items()
        if k not in ("grid_lookup", "station_dir")
    }
    config_sha = hashlib.sha256(
        json.dumps(
            {**portable, "defo_db_version": defo_db_version}, sort_keys=True
        ).encode()
    ).hexdigest()
    provenance = Provenance(
        grid_version=cfg.version,
        grid_type=cfg.grid_type,
        reference_frame=cfg.reference_frame,
        snapshot_id=cfg.snapshot_id,
        lookup_sha256=_sha256_file(Path(cfg.grid_lookup)),
        buffer_meters=cfg.buffer_meters,
        n_nodes=len(nodes),
        n_missing=n_missing,
        n_excluded=n_excluded,
        n_reinterpolated=n_reinterpolated,
        defo_db_version=defo_db_version,
        field_sha256=_sha256_arrays(*(fields[c] for p in VALUE_COLUMNS for c in p)),
        config_sha256=config_sha,
    )
    logger.info(
        "Sampled %d UNR grid nodes (%s, %s, buffer %.0f m) onto a %s grid; digest %s",
        len(nodes),
        cfg.grid_type,
        cfg.reference_frame,
        cfg.buffer_meters,
        fields["ve"].shape,
        provenance.digest[:12],
    )
    return GnssField(
        ve=fields["ve"],
        vn=fields["vn"],
        vu=fields["vu"],
        sigma_ve=fields["sigma_ve"],
        sigma_vn=fields["sigma_vn"],
        sigma_vu=fields["sigma_vu"],
        nodes=nodes,
        provenance=provenance,
    )


def project_field_to_los(
    gnss: GnssField,
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    dt_years: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Project the gridded E/N/U field onto the LOS, per pixel, scaled by `dt_years`.

    Returns ``(los, sigma_los)`` in the field's units (mm/yr x years = mm).
    Pixels without a look (``e == n == 0`` or NaN LOS) are NaN. Sigmas are
    combined assuming uncorrelated components (the UNR tenv8 files carry no
    correlations).
    """
    e, n, u = (np.asarray(c, dtype=np.float64) for c in (los_east, los_north, los_up))
    if e.shape != gnss.grid_shape:
        msg = f"LOS rasters {e.shape} do not match the GNSS field {gnss.grid_shape}"
        raise ValueError(msg)
    no_look = (
        ((e == 0) & (n == 0)) | ~np.isfinite(e) | ~np.isfinite(n) | ~np.isfinite(u)
    )
    los = (e * gnss.ve + n * gnss.vn + u * gnss.vu) * dt_years
    sigma = np.sqrt(
        (e * gnss.sigma_ve) ** 2 + (n * gnss.sigma_vn) ** 2 + (u * gnss.sigma_vu) ** 2
    ) * abs(dt_years)
    los = np.where(no_look, np.nan, los)
    sigma = np.where(no_look, np.nan, sigma)
    return los.astype(np.float32), sigma.astype(np.float32)
