# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""GNSS E/N/U for VLM, consistent with the DISP-CAL products it uses (plan T55).

VLM decomposes calibrated DISP into vertical (and east) motion and takes the
north component from GNSS (PRD §3.2, D15). That GNSS must be the same
realisation DISP-CAL was calibrated to, or VLM mixes two reference frames.
`load_gnss_for_vlm`:

1. reads the GNSS provenance record of every DISP-CAL product
   (``/metadata/gnss_provenance``, a JSON string; plan T37.3);
2. takes the GNSS field from the products' ``gnss_ve``/``gnss_vn`` layers
   when they carry them (plan T61), otherwise samples the grid snapshot
   with `sample_gnss_enu`, the single sampling path;
3. verifies that every product and the field agree on the realisation
   (`verify_provenance`) and raises `ProvenanceMismatchError` if not.

Products written before T37.3 carry no record: `ProvenanceMissingError` is raised
unless ``allow_missing_provenance=True``, which logs a warning instead.
`pair_displacement` scales the velocities to a DISP pair (T55.3).
"""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping, Sequence
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np

from .sampling import GnssField, GnssGridConfig, Provenance, sample_gnss_enu

logger = logging.getLogger(__name__)

__all__ = [
    "PROVENANCE_KEYS",
    "ProvenanceMismatchError",
    "ProvenanceMissingError",
    "load_gnss_for_vlm",
    "pair_displacement",
    "read_cal_provenance",
    "verify_provenance",
]

# What defines the realisation. Node counts and the field hash depend on the
# target grid (VLM may run on another grid), so they are not compared.
PROVENANCE_KEYS: tuple[str, ...] = (
    "grid_version",
    "grid_type",
    "reference_frame",
    "snapshot_id",
    "lookup_sha256",
    "buffer_meters",
    "defo_db_version",
    "config_sha256",
)


class ProvenanceMismatchError(ValueError):
    """The GNSS realisation differs between DISP-CAL and VLM."""


class ProvenanceMissingError(ValueError):
    """A DISP-CAL product carries no GNSS provenance record."""


def _as_mapping(p: Provenance | Mapping[str, Any]) -> Mapping[str, Any]:
    return p.as_dict() if isinstance(p, Provenance) else p


def verify_provenance(
    expected: Provenance | Mapping[str, Any],
    actual: Provenance | Mapping[str, Any],
    *,
    source: str = "",
) -> None:
    """Raise `ProvenanceMismatchError` unless every `PROVENANCE_KEYS` entry agrees.

    The message names each differing key with both values, e.g. a grid
    version or a defo-area database version that moved between the DISP-CAL
    run and the VLM run.
    """
    e, a = _as_mapping(expected), _as_mapping(actual)
    diffs = [
        f"{k}: DISP-CAL {e.get(k)!r} vs VLM {a.get(k)!r}"
        for k in PROVENANCE_KEYS
        if e.get(k) != a.get(k)
    ]
    if diffs:
        where = f" ({source})" if source else ""
        msg = "GNSS realisation differs from DISP-CAL" + where + ": " + "; ".join(diffs)
        raise ProvenanceMismatchError(msg)


def _read_string(obj: Any) -> str:
    value = obj[()] if hasattr(obj, "shape") else obj
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        value = value.item()
        return value.decode() if isinstance(value, bytes) else str(value)
    return str(value)


def read_cal_provenance(path: Path | str) -> dict[str, Any] | None:
    """GNSS provenance record of a DISP-CAL product, or None if it has none.

    Read from ``/metadata/gnss_provenance`` (dataset or attribute, a JSON
    string).
    """
    import h5py

    with h5py.File(path, "r") as f:
        meta = f.get("metadata")
        if meta is None:
            return None
        if "gnss_provenance" in meta:
            text = _read_string(meta["gnss_provenance"])
        elif "gnss_provenance" in meta.attrs:
            text = _read_string(meta.attrs["gnss_provenance"])
        else:
            return None
    return json.loads(text)


def _read_cal_gnss_layers(path: Path | str) -> dict[str, np.ndarray] | None:
    """Read the ``gnss_ve``/``gnss_vn`` (+ sigma, optional up) layers (plan T61)."""
    import h5py

    names = (
        "gnss_ve",
        "gnss_vn",
        "gnss_vu",
        "gnss_ve_std",
        "gnss_vn_std",
        "gnss_vu_std",
    )
    with h5py.File(path, "r") as f:
        if "gnss_ve" not in f or "gnss_vn" not in f:
            return None
        return {n: np.squeeze(f[n][()]).astype(np.float64) for n in names if n in f}


def load_gnss_for_vlm(
    cal_products: Sequence[Path | str],
    grid_cfg: GnssGridConfig,
    grid: Any,
    *,
    exclude: Any | None = None,
    defo_db_version: str | None = None,
    allow_missing_provenance: bool = False,
) -> GnssField:
    """GNSS E/N/U for VLM, verified against the DISP-CAL products.

    Parameters
    ----------
    cal_products : sequence of path
        The DISP-CAL products VLM uses (ascending and descending, all pairs).
    grid_cfg, grid, exclude, defo_db_version
        As `sample_gnss_enu`; used when the products carry no GNSS layers.
    allow_missing_provenance : bool
        Accept products without a provenance record (written before plan
        T37.3) with a warning instead of `ProvenanceMissingError`.

    Raises
    ------
    ProvenanceMissingError
        A product has no record and `allow_missing_provenance` is False.
    ProvenanceMismatchError
        Two products, or a product and the sampled field, differ in any of
        `PROVENANCE_KEYS`.

    """
    if not cal_products:
        msg = "no DISP-CAL products given"
        raise ValueError(msg)
    records: list[tuple[str, dict[str, Any]]] = []
    for p in cal_products:
        rec = read_cal_provenance(p)
        if rec is None:
            if not allow_missing_provenance:
                msg = (
                    f"{Path(p).name} has no /metadata/gnss_provenance; it predates"
                    " the provenance record (plan T37.3). Pass"
                    " allow_missing_provenance=True to accept it unverified."
                )
                raise ProvenanceMissingError(msg)
            logger.warning("%s: no GNSS provenance record; not verified", Path(p).name)
            continue
        records.append((Path(p).name, rec))
    for name, rec in records[1:]:
        verify_provenance(records[0][1], rec, source=f"{records[0][0]} vs {name}")

    layers = _read_cal_gnss_layers(cal_products[0])
    if layers is not None and records:
        field = _field_from_layers(layers, records[0][1])
        logger.info(
            "GNSS E/N from the DISP-CAL layers of %s", Path(cal_products[0]).name
        )
    else:
        field = sample_gnss_enu(grid_cfg, grid, exclude, defo_db_version)
        logger.info(
            "GNSS E/N sampled from the grid snapshot (digest %s)",
            field.provenance.digest[:12],
        )
    for name, rec in records:
        verify_provenance(rec, field.provenance, source=name)
    return field


def _field_from_layers(
    layers: Mapping[str, np.ndarray], record: Mapping[str, Any]
) -> GnssField:
    import pandas as pd

    nan = np.full(layers["gnss_ve"].shape, np.nan)
    keys = Provenance.__dataclass_fields__
    prov = Provenance(**{k: record[k] for k in keys})
    return GnssField(
        ve=layers["gnss_ve"],
        vn=layers["gnss_vn"],
        vu=layers.get("gnss_vu", nan),
        sigma_ve=layers.get("gnss_ve_std", nan),
        sigma_vn=layers.get("gnss_vn_std", nan),
        sigma_vu=layers.get("gnss_vu_std", nan),
        nodes=pd.DataFrame(),
        provenance=prov,
    )


def _decimal_year(d: date | datetime) -> float:
    d0 = datetime(d.year, 1, 1)
    d1 = datetime(d.year + 1, 1, 1)
    dt = datetime(d.year, d.month, d.day) if not isinstance(d, datetime) else d
    return d.year + (dt - d0).total_seconds() / (d1 - d0).total_seconds()


def pair_displacement(
    field: GnssField, reference: date | datetime, secondary: date | datetime
) -> dict[str, np.ndarray]:
    """E/N/U displacement (and sigma) of the pair: velocity x interval (T55.3).

    Returns ``{"de", "dn", "du", "sigma_de", "sigma_dn", "sigma_du"}`` in the
    field's units times years (mm for the UNR grid). The constant grid holds
    secular rates only, so this is the secular displacement over the pair.
    """
    dt = _decimal_year(secondary) - _decimal_year(reference)
    return {
        "de": field.ve * dt,
        "dn": field.vn * dt,
        "du": field.vu * dt,
        "sigma_de": field.sigma_ve * abs(dt),
        "sigma_dn": field.sigma_vn * abs(dt),
        "sigma_du": field.sigma_vu * abs(dt),
    }
