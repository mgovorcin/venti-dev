# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""sample_gnss_enu: the single GNSS sampling path (plan T28).

A synthetic UNR grid: a planar E/N/U velocity field, nodes on a 10 km lattice
inside and around a 40 x 40 km UTM frame, constant-grid tenv8 files written
in the layout `venti.gnss.unr` uses. Truth is the plane, so the gridded field
can be checked exactly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from pyproj import Transformer

from venti.gnss.sampling import GnssGridConfig, project_field_to_los, sample_gnss_enu

UTM = 32615  # Houston
X0, Y0 = 250_000.0, 3_300_000.0  # frame lower-left, metres
SIZE = 40_000.0
NX = NY = 40  # 1 km posting


def plane(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """mm/yr fields linear in UTM coordinates (relative to the frame centre)."""
    dx = (x - (X0 + SIZE / 2)) / 1e5
    dy = (y - (Y0 + SIZE / 2)) / 1e5
    return -12.0 + 3.0 * dx, -2.0 - 1.0 * dy, -3.0 + 2.0 * dx + 1.5 * dy


@pytest.fixture
def grid() -> tuple[np.ndarray, np.ndarray]:
    x = X0 + (np.arange(NX) + 0.5) * (SIZE / NX)
    y = Y0 + SIZE - (np.arange(NY) + 0.5) * (SIZE / NY)
    return x, y


@pytest.fixture
def snapshot(tmp_path: Path) -> tuple[Path, Path, np.ndarray]:
    """Lookup + constant tenv8 files for nodes on a 10 km lattice, 30 km beyond the frame."""
    to_ll = Transformer.from_crs(f"EPSG:{UTM}", "EPSG:4326", always_xy=True)
    xs = np.arange(X0 - 30_000, X0 + SIZE + 30_001, 10_000.0)
    ys = np.arange(Y0 - 30_000, Y0 + SIZE + 30_001, 10_000.0)
    xx, yy = np.meshgrid(xs, ys)
    xx, yy = xx.ravel(), yy.ravel()
    ids = np.arange(1, xx.size + 1)
    lon, lat = to_ll.transform(xx, yy)
    lookup = tmp_path / "grid_latlon_lookup.txt"
    np.savetxt(lookup, np.c_[ids, lon, lat], fmt=["%06d", "%.8f", "%.8f"])
    station_dir = tmp_path / "gnss"
    station_dir.mkdir()
    ve, vn, vu = plane(xx, yy)
    t = 2015 + np.arange(0, 8) * 0.5
    for i, node in enumerate(ids):
        # year east north up sigma_e sigma_n sigma_u (+ a dummy column)
        rows = np.c_[
            t,
            ve[i] * (t - 2015),
            vn[i] * (t - 2015),
            vu[i] * (t - 2015),
            np.full_like(t, 0.4),
            np.full_like(t, 0.5),
            np.full_like(t, 1.2),
            np.zeros_like(t),
        ]
        np.savetxt(station_dir / f"{node:06d}_IGS20_constant.tenv8", rows, fmt="%.6f")
    return lookup, station_dir, np.c_[xx, yy]


def interior_ids(xy: np.ndarray) -> np.ndarray:
    """1-based ids of lattice nodes at least 5 km inside the frame."""
    inside = (
        (xy[:, 0] > X0 + 5_000)
        & (xy[:, 0] < X0 + SIZE - 5_000)
        & (xy[:, 1] > Y0 + 5_000)
        & (xy[:, 1] < Y0 + SIZE - 5_000)
    )
    return np.flatnonzero(inside) + 1


def make_cfg(snapshot, **kw) -> GnssGridConfig:
    lookup, station_dir, _ = snapshot
    return GnssGridConfig(
        grid_lookup=lookup, station_dir=station_dir, utm_epsg=UTM, **kw
    )


def test_recovers_the_plane_and_sigmas(snapshot, grid):
    gnss = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    gx, gy = np.meshgrid(*grid)
    ve, vn, vu = plane(gx, gy)
    assert gnss.grid_shape == (NY, NX)
    np.testing.assert_allclose(gnss.ve, ve, atol=0.05)
    np.testing.assert_allclose(gnss.vn, vn, atol=0.05)
    np.testing.assert_allclose(gnss.vu, vu, atol=0.05)
    np.testing.assert_allclose(gnss.sigma_vu, 1.2, atol=0.02)
    assert (gnss.sigma_ve >= 0).all()


def test_buffer_adds_nodes_outside_the_frame(snapshot, grid):
    inside_only = sample_gnss_enu(make_cfg(snapshot), grid)
    buffered = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    assert buffered.provenance.n_nodes > inside_only.provenance.n_nodes
    assert inside_only.provenance.buffer_meters == 0.0
    # without the buffer the field is an extrapolation at the frame edge:
    # the buffered one is closer to the truth there
    gx, gy = np.meshgrid(*grid)
    ve, _, _ = plane(gx, gy)
    edge = np.s_[0, :]
    assert (
        np.abs(buffered.ve[edge] - ve[edge]).max()
        <= np.abs(inside_only.ve[edge] - ve[edge]).max() + 1e-6
    )


def test_variable_grid_is_gated(snapshot, grid):
    with pytest.raises(ValueError, match="reprocessing=True"):
        sample_gnss_enu(make_cfg(snapshot, grid_type="variable"), grid)


def test_missing_station_files_are_skipped_and_counted(snapshot, grid):
    _lookup, station_dir, xy = snapshot
    before = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    for node in interior_ids(xy)[:3]:
        (station_dir / f"{node:06d}_IGS20_constant.tenv8").unlink()
    gnss = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    assert gnss.provenance.n_missing == 3
    assert gnss.provenance.n_nodes == before.provenance.n_nodes - 3
    assert before.provenance.n_missing == 0


def test_exclusion_drops_or_reinterpolates(snapshot, grid):
    shapely = pytest.importorskip("shapely")
    _lookup, station_dir, xy = snapshot
    # corrupt the nodes inside a box (a "subsidence bowl" in the grid) by -20 mm/yr up
    to_ll = Transformer.from_crs(f"EPSG:{UTM}", "EPSG:4326", always_xy=True)
    box_utm = (X0 + 10_000, Y0 + 10_000, X0 + 25_000, Y0 + 25_000)
    lon_min, lat_min = to_ll.transform(box_utm[0], box_utm[1])
    lon_max, lat_max = to_ll.transform(box_utm[2], box_utm[3])
    area = shapely.box(lon_min - 1e-3, lat_min - 1e-3, lon_max + 1e-3, lat_max + 1e-3)
    lon_all, lat_all = to_ll.transform(xy[:, 0], xy[:, 1])
    inside = np.asarray(shapely.contains_xy(area, lon_all, lat_all), dtype=bool)
    for i in np.flatnonzero(inside):
        f = station_dir / f"{i + 1:06d}_IGS20_constant.tenv8"
        rows = np.loadtxt(f)
        rows[:, 3] -= 20.0 * (rows[:, 0] - 2015)  # extra -20 mm/yr up
        np.savetxt(f, rows, fmt="%.6f")
    gx, gy = np.meshgrid(*grid)
    _, _, vu_truth = plane(gx, gy)

    polluted = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    assert np.nanmin(polluted.vu - vu_truth) < -5  # the bowl leaked into the field

    dropped = sample_gnss_enu(
        make_cfg(snapshot, buffer_meters=30_000, exclude_defo_nodes=True),
        grid,
        exclude=area,
        defo_db_version="v1",
    )
    assert dropped.provenance.n_excluded == int(inside.sum()) > 0
    assert dropped.provenance.n_reinterpolated == 0
    assert dropped.provenance.defo_db_version == "v1"
    np.testing.assert_allclose(dropped.vu, vu_truth, atol=0.3)

    reinterp = sample_gnss_enu(
        make_cfg(
            snapshot,
            buffer_meters=30_000,
            exclude_defo_nodes=True,
            reinterpolate_excluded=True,
        ),
        grid,
        exclude=area,
    )
    assert reinterp.provenance.n_reinterpolated == int(inside.sum())
    assert reinterp.provenance.n_nodes == polluted.provenance.n_nodes  # nodes kept
    np.testing.assert_allclose(reinterp.vu, vu_truth, atol=1.0)
    assert (
        reinterp.provenance.digest
        != dropped.provenance.digest
        != polluted.provenance.digest
    )


def test_provenance_tracks_config_and_data(snapshot, grid):
    a = sample_gnss_enu(make_cfg(snapshot, snapshot_id="unr-0.3-2026-09"), grid)
    b = sample_gnss_enu(make_cfg(snapshot, snapshot_id="unr-0.3-2026-09"), grid)
    assert a.provenance == b.provenance  # deterministic
    assert a.provenance.as_dict()["digest"] == a.provenance.digest
    assert len(a.provenance.lookup_sha256) == 64
    c = sample_gnss_enu(make_cfg(snapshot, snapshot_id="unr-0.3-2027-03"), grid)
    assert c.provenance.field_sha256 == a.provenance.field_sha256  # same data ...
    assert c.provenance.digest != a.provenance.digest  # ... different realisation
    # a changed node changes the field hash
    _lookup, station_dir, xy = snapshot
    f = station_dir / f"{interior_ids(xy)[0]:06d}_IGS20_constant.tenv8"
    rows = np.loadtxt(f)
    rows[:, 1] += 1.0 * (rows[:, 0] - 2015)
    np.savetxt(f, rows, fmt="%.6f")
    d = sample_gnss_enu(make_cfg(snapshot, snapshot_id="unr-0.3-2026-09"), grid)
    assert d.provenance.field_sha256 != a.provenance.field_sha256


def test_project_field_to_los(snapshot, grid):
    gnss = sample_gnss_enu(make_cfg(snapshot, buffer_meters=30_000), grid)
    e = np.full(gnss.grid_shape, -0.6)
    n = np.full(gnss.grid_shape, 0.1)
    u = np.full(gnss.grid_shape, 0.79)
    u[0, 0] = 1.0  # no-look pixel marked (0, 0, 1)
    e[0, 0] = n[0, 0] = 0.0
    los, sigma = project_field_to_los(gnss, e, n, u, dt_years=0.5)
    expected = (e * gnss.ve + n * gnss.vn + u * gnss.vu) * 0.5
    np.testing.assert_allclose(los[1:, 1:], expected[1:, 1:], rtol=1e-5)
    assert np.isnan(los[0, 0])
    assert np.isnan(sigma[0, 0])
    assert los.dtype == np.float32
    assert (sigma[1:, 1:] > 0).all()
    with pytest.raises(ValueError, match="do not match"):
        project_field_to_los(gnss, e[:5], n[:5], u[:5])


def test_no_nodes_in_bounds_raises(snapshot):
    far = (np.array([900_000.0, 901_000.0]), np.array([4_900_000.0, 4_899_000.0]))
    with pytest.raises(ValueError, match="No UNR grid nodes"):
        sample_gnss_enu(make_cfg(snapshot), far)


def test_exclusion_with_unr_0_360_longitudes(snapshot, grid):
    """Regression: UNR lookups store 0-360 longitudes; the exclusion compared
    them with -180..180 polygons and never matched (R-G5 was a no-op)."""
    shapely = pytest.importorskip("shapely")
    lookup, _station_dir, xy = snapshot
    ids, lon, lat = np.loadtxt(lookup, unpack=True)
    np.savetxt(
        lookup, np.c_[ids, np.mod(lon, 360.0), lat], fmt=["%06d", "%.8f", "%.8f"]
    )
    to_ll = Transformer.from_crs(f"EPSG:{UTM}", "EPSG:4326", always_xy=True)
    lo0, la0 = to_ll.transform(X0 + 10_000, Y0 + 10_000)
    lo1, la1 = to_ll.transform(X0 + 25_000, Y0 + 25_000)
    area = shapely.box(lo0 - 1e-3, la0 - 1e-3, lo1 + 1e-3, la1 + 1e-3)
    lon_all, lat_all = to_ll.transform(xy[:, 0], xy[:, 1])
    n_inside = int(shapely.contains_xy(area, lon_all, lat_all).sum())
    out = sample_gnss_enu(
        make_cfg(snapshot, buffer_meters=30_000, exclude_defo_nodes=True),
        grid,
        exclude=area,
    )
    assert out.provenance.n_excluded == n_inside > 0
    assert out.nodes["lon"].between(-180, 180).all()
