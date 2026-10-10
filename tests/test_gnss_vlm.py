# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""GNSS for VLM with the DISP-CAL provenance check (plan T55)."""

from __future__ import annotations

import json
import shutil
from datetime import date
from pathlib import Path

import h5py
import numpy as np
import pytest

from venti.gnss.sampling import sample_gnss_enu
from venti.gnss.vlm import (
    ProvenanceMismatchError,
    ProvenanceMissingError,
    load_gnss_for_vlm,
    pair_displacement,
    read_cal_provenance,
    verify_provenance,
)

from .test_gnss_sampling import make_cfg


def _cal_product(
    path: Path, provenance: dict | None, layers: dict | None = None
) -> Path:
    """A minimal DISP-CAL file: /metadata/gnss_provenance and optional layers."""
    with h5py.File(path, "w") as f:
        meta = f.create_group("metadata")
        meta.create_dataset("cal_disp_software_version", data="0.5.0")
        if provenance is not None:
            meta.create_dataset("gnss_provenance", data=json.dumps(provenance))
        for name, arr in (layers or {}).items():
            f[name] = arr
    return path


def test_matching_provenance_passes_and_returns_the_sampled_field(
    tmp_path, snapshot, grid
):
    cfg = make_cfg(snapshot, buffer_meters=30_000, snapshot_id="unr-0.3")
    ref = sample_gnss_enu(cfg, grid)  # what DISP-CAL recorded
    prods = [
        _cal_product(tmp_path / f"cal{i}.nc", ref.provenance.as_dict())
        for i in range(2)
    ]
    field = load_gnss_for_vlm(prods, cfg, grid)
    np.testing.assert_array_equal(field.vn, ref.vn)
    assert read_cal_provenance(prods[0])["digest"] == ref.provenance.digest


def test_portable_hash_same_snapshot_elsewhere(tmp_path, snapshot, grid):
    """The same snapshot at another path is the same realisation."""
    lookup, station_dir, xy = snapshot
    moved = tmp_path / "elsewhere"
    shutil.copytree(station_dir, moved / "gnss")
    shutil.copy(lookup, moved / lookup.name)
    a = sample_gnss_enu(make_cfg(snapshot), grid)
    b = sample_gnss_enu(make_cfg((moved / lookup.name, moved / "gnss", xy)), grid)
    assert a.provenance.config_sha256 == b.provenance.config_sha256
    verify_provenance(a.provenance, b.provenance)


@pytest.mark.parametrize(
    ("change", "key"),
    [
        ({"snapshot_id": "unr-0.3-other"}, "snapshot_id"),
        ({"buffer_meters": 50_000}, "buffer_meters"),
    ],
)
def test_mismatched_realisation_raises(tmp_path, snapshot, grid, change, key):
    cal_cfg = make_cfg(snapshot, buffer_meters=30_000, snapshot_id="unr-0.3")
    rec = sample_gnss_enu(cal_cfg, grid).provenance.as_dict()
    prod = _cal_product(tmp_path / "cal.nc", rec)
    vlm_cfg = make_cfg(
        snapshot, **{"buffer_meters": 30_000, "snapshot_id": "unr-0.3", **change}
    )
    with pytest.raises(ProvenanceMismatchError, match=key):
        load_gnss_for_vlm([prod], vlm_cfg, grid)


def test_grid_and_defo_db_version_mismatch_messages():
    a = {"grid_version": "0.3", "defo_db_version": "v1", "lookup_sha256": "x"}
    with pytest.raises(
        ProvenanceMismatchError, match=r"grid_version: DISP-CAL '0.3' vs VLM '0.4'"
    ):
        verify_provenance(a, {**a, "grid_version": "0.4"}, source="cal.nc")
    with pytest.raises(ProvenanceMismatchError, match=r"defo_db_version.*'v1'.*'v2'"):
        verify_provenance(a, {**a, "defo_db_version": "v2"})
    verify_provenance(a, dict(a))


def test_products_must_agree_with_each_other(tmp_path, snapshot, grid):
    cfg = make_cfg(snapshot, snapshot_id="unr-0.3")
    rec = sample_gnss_enu(cfg, grid).provenance.as_dict()
    asc = _cal_product(tmp_path / "asc.nc", rec)
    desc = _cal_product(tmp_path / "desc.nc", {**rec, "grid_type": "variable"})
    with pytest.raises(ProvenanceMismatchError, match=r"asc\.nc vs desc\.nc"):
        load_gnss_for_vlm([asc, desc], cfg, grid)


def test_missing_provenance_fails_unless_allowed(tmp_path, snapshot, grid, caplog):
    cfg = make_cfg(snapshot)
    old = _cal_product(tmp_path / "gamma.nc", None)
    with pytest.raises(ProvenanceMissingError, match=r"T37\.3"):
        load_gnss_for_vlm([old], cfg, grid)
    field = load_gnss_for_vlm([old], cfg, grid, allow_missing_provenance=True)
    assert np.isfinite(field.vn).all()
    assert "not verified" in caplog.text


def test_layers_from_the_product_are_used(tmp_path, snapshot, grid):
    """T61 layers present: no resampling, the product's E/N are returned."""
    cfg = make_cfg(snapshot, snapshot_id="unr-0.3")
    rec = sample_gnss_enu(cfg, grid).provenance.as_dict()
    shape = (len(grid[1]), len(grid[0]))
    layers = {"gnss_ve": np.full(shape, -11.0), "gnss_vn": np.full(shape, -2.5)}
    prod = _cal_product(tmp_path / "cal.nc", rec, layers)
    field = load_gnss_for_vlm([prod], cfg, grid)
    assert np.all(field.vn == -2.5)
    assert np.isnan(field.vu).all()
    assert field.provenance.digest == rec["digest"]


def test_pair_displacement_scales_by_the_interval(snapshot, grid):
    field = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    d = pair_displacement(field, date(2022, 1, 1), date(2022, 7, 2))
    dt = d["dn"] / field.vn
    np.testing.assert_allclose(dt, 182 / 365, rtol=1e-3)
    np.testing.assert_allclose(d["sigma_dn"], field.sigma_vn * dt)
    back = pair_displacement(field, date(2022, 7, 2), date(2022, 1, 1))
    np.testing.assert_allclose(back["dn"], -d["dn"])
    np.testing.assert_allclose(back["sigma_dn"], d["sigma_dn"])
