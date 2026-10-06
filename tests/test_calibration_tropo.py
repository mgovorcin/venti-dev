# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tropo modes and the cal-disp parity (PRD R-T1; plan T34).

The parity tests import ``cal_disp`` and are skipped when it is not on the
path (``PYTHONPATH=<cal-disp>/src pytest tests/test_calibration_tropo.py``).
The golden-pair independent check needs ``CAL_DISP_GOLDEN_DIR``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator

from venti.calibration import tropo
from venti.calibration.two_pass import calibrate_pair
from venti.frames import load_frame_table
from venti.workflow.config import (
    AlgorithmParameters,
    CalibrationOptions,
    TropoOptions,
)

UTM = 32615
NY = NX = 256
PIXEL_M = 300.0
X0, Y0 = 250_000.0, 3_300_000.0 + NY * PIXEL_M  # upper-left corner


@pytest.fixture(scope="module")
def dem() -> xr.DataArray:
    """256 x 256 DEM in UTM 15N: a 1.2 km ridge on a gentle slope, with NaN holes."""
    rng = np.random.default_rng(1)
    x = X0 + (np.arange(NX) + 0.5) * PIXEL_M
    y = Y0 - (np.arange(NY) + 0.5) * PIXEL_M
    xx, yy = np.meshgrid(x, y)
    h = (
        50.0
        + 0.002 * (xx - X0)
        + 1200.0 * np.exp(-(((xx - X0 - 40e3) / 12e3) ** 2))
        + 5.0 * rng.standard_normal((NY, NX))
    )
    h[10:20, 200:230] = np.nan  # a void, as in the STATIC DEM over water
    da = xr.DataArray(h.astype(np.float32), coords={"y": y, "x": x}, dims=("y", "x"))
    return da.rio.write_crs(UTM)


@pytest.fixture(scope="module")
def cube() -> xr.DataArray:
    """Zenith delay cube in EPSG:4326 covering the DEM: exponential in height
    plus a lateral gradient and a lateral 'turbulent' bump."""
    heights = np.array([0, 250, 500, 1000, 1500, 2000, 3000, 4500, 6000, 9000.0])
    lat = np.linspace(29.0, 31.5, 26)
    lon = np.linspace(-96.5, -94.0, 26)
    hh, la, lo = np.meshgrid(heights, lat, lon, indexing="ij")
    delay = (
        2.4 * np.exp(-hh / 8000.0)
        + 0.02 * (lo + 95.25)
        + 0.015 * np.exp(-(((la - 30.1) / 0.3) ** 2 + ((lo + 95.3) / 0.3) ** 2))
    )
    da = xr.DataArray(
        delay.astype(np.float32),
        coords={"height": heights, "latitude": lat, "longitude": lon},
        dims=("height", "latitude", "longitude"),
        name="zenith_total_delay",
    )
    return da


@pytest.fixture(scope="module")
def los_up(dem) -> xr.DataArray:
    return xr.full_like(dem, 0.8, dtype=np.float32).rio.write_crs(UTM)


# --------------------------------------------------------------------------
# parity with cal-disp


@dataclass
class _FakeProduct:
    date: datetime
    filename: str
    da: xr.DataArray

    def get_total_delay(self, **_kw):
        return self.da


class TestParityWithCalDisp:
    """T34.1: the ported numerics are bit-identical to cal_disp.product._tropo."""

    @pytest.fixture(autouse=True)
    def _cal_disp(self):
        self.cd = pytest.importorskip("cal_disp.product._tropo")

    def test_interpolate_to_dem_surface_is_bit_identical(self, cube, dem):
        ours = tropo.interpolate_to_dem_surface(cube, dem)
        theirs = self.cd.interpolate_to_dem_surface(cube, dem)
        assert ours.dtype == theirs.dtype == np.float32
        np.testing.assert_array_equal(ours.values, theirs.values)
        assert ours.attrs == theirs.attrs
        # and the block size does not matter (the 512-row loop is memory only)
        blocked = tropo.interpolate_to_dem_surface(cube, dem, block_rows=37)
        np.testing.assert_array_equal(ours.values, blocked.values)

    def test_compute_los_correction_is_bit_identical(self, cube, dem, los_up):
        ztd = tropo.interpolate_to_dem_surface(cube, dem)
        ours = tropo.compute_los_correction(ztd, los_up)
        theirs = self.cd.compute_los_correction(ztd, los_up)
        np.testing.assert_array_equal(ours.values, theirs.values)
        ref = ours * 0.9
        ours2 = tropo.compute_los_correction(ztd, los_up, reference_correction=ref)
        theirs2 = self.cd.compute_los_correction(ztd, los_up, reference_correction=ref)
        np.testing.assert_array_equal(ours2.values, theirs2.values)
        assert ours2.dtype == np.float32

    def test_interpolate_in_time_is_bit_identical(self, cube):
        t0, t1, t = datetime(2022, 1, 11, 0), datetime(2022, 1, 11, 6), datetime(
            2022, 1, 11, 0, 17
        )
        late = cube * 1.03
        ours = tropo.interpolate_in_time(
            cube, late, t0, t1, t, early_name="a.nc", late_name="b.nc"
        )
        theirs = self.cd.interpolate_in_time(
            _FakeProduct(t0, "a.nc", cube), _FakeProduct(t1, "b.nc", late), t
        )
        np.testing.assert_array_equal(ours.values, theirs.values)
        assert ours.attrs == theirs.attrs


# --------------------------------------------------------------------------
# the numerics on their own


class TestNumerics:
    def test_surface_delay_matches_the_analytic_profile(self, cube, dem):
        ztd = tropo.interpolate_to_dem_surface(cube, dem)
        h = dem.values.astype(float)
        ok = np.isfinite(h)
        expected = 2.4 * np.exp(-h / 8000.0)
        # lateral terms are < 0.04 m; the height profile is linear-interpolated
        # between levels 250 m apart near the ground
        err = (ztd.values - expected)[ok]
        assert np.nanmax(np.abs(err)) < 0.05
        assert np.isnan(ztd.values[~ok]).all()
        assert ztd.rio.crs == dem.rio.crs

    def test_time_interpolation_rules(self, cube):
        t0, t1 = datetime(2022, 1, 11, 0), datetime(2022, 1, 11, 6)
        mid = tropo.interpolate_in_time(cube, cube * 3, t0, t1, datetime(2022, 1, 11, 3))
        np.testing.assert_allclose(mid.values, 2 * cube.values, rtol=1e-6)
        assert mid.attrs["interpolation_weight"] == 0.5
        with pytest.raises(ValueError, match="must be before"):
            tropo.interpolate_in_time(cube, cube, t1, t0, t0)
        with pytest.raises(ValueError, match="between"):
            tropo.interpolate_in_time(cube, cube, t0, t1, datetime(2022, 1, 12))

    def test_los_sign_and_pair_difference(self, cube, dem, los_up):
        ztd = tropo.interpolate_to_dem_surface(cube, dem)
        los = tropo.compute_los_correction(ztd, los_up)
        np.testing.assert_allclose(los.values, -ztd.values / 0.8, rtol=1e-6)
        corr = tropo.pair_correction(los.values, 1.1 * los.values)
        assert corr.dtype == np.float32
        np.testing.assert_allclose(corr, 0.1 * los.values, rtol=1e-5)
        with pytest.raises(ValueError, match="differ"):
            tropo.pair_correction(los.values, los.values[:10])
        with pytest.raises(ValueError, match="CRS"):
            tropo.interpolate_to_dem_surface(cube, xr.DataArray(dem.values))


# --------------------------------------------------------------------------
# stratified fit


def _grid_xy():
    x = X0 + (np.arange(NX) + 0.5) * PIXEL_M
    y = Y0 - (np.arange(NY) + 0.5) * PIXEL_M
    return x, y


class TestStratified:
    def test_recovers_height_terms_and_drops_the_turbulent_part(self, dem):
        rng = np.random.default_rng(2)
        h = dem.values.astype(float)
        x, y = _grid_xy()
        xx, yy = np.meshgrid(x, y)
        strat = 0.001 - 0.006 * (h / 1e3) + 0.0008 * (h / 1e3) ** 2  # -6 mm/km
        plane = 1e-8 * (xx - X0) - 2e-8 * (yy - Y0 + NY * PIXEL_M)
        turb = 0.004 * np.sin(xx / 1500.0) * np.cos(yy / 1100.0)  # ~9 km lateral noise
        delay = strat + plane + turb + 0.0002 * rng.standard_normal(h.shape)
        model = tropo.stratified_tropo(delay, h, x=x, y=y)
        assert model.gradient_per_km == pytest.approx(-0.006, abs=0.0005)
        assert model.n_fit > 1000
        fitted = tropo.stratified_delay(delay, h, x=x, y=y)
        ok = np.isfinite(h)
        np.testing.assert_allclose(fitted[ok], (strat + plane)[ok], atol=0.0007)
        assert np.isnan(fitted[~ok]).all()
        # the turbulent part is what is left over, not what is applied
        resid = (delay - fitted)[ok]
        assert np.corrcoef(resid, turb[ok])[0, 1] > 0.95

    def test_fit_is_linear_and_parametrisation_free(self, dem):
        h = dem.values.astype(float)
        x, y = _grid_xy()
        rng = np.random.default_rng(3)
        ref = -0.005 * h / 1e3 + 0.001 * rng.standard_normal(h.shape)
        sec = -0.008 * h / 1e3 + 0.001 * rng.standard_normal(h.shape)
        pair = tropo.stratified_delay(sec - ref, h, x=x, y=y)
        per_date = tropo.stratified_delay(sec, h, x=x, y=y) - tropo.stratified_delay(
            ref, h, x=x, y=y
        )
        ok = np.isfinite(h)
        np.testing.assert_allclose(pair[ok], per_date[ok], atol=1e-9)
        # pixel indices instead of map coordinates span the same column space
        idx = tropo.stratified_delay(sec - ref, h)
        np.testing.assert_allclose(pair[ok], idx[ok], atol=1e-9)

    def test_valid_mask_and_errors(self, dem):
        h = dem.values.astype(float)
        delay = -0.005 * h / 1e3
        valid = np.ones(h.shape, bool)
        valid[:, 128:] = False
        model = tropo.stratified_tropo(delay, h, valid)
        assert model.gradient_per_km == pytest.approx(-0.005, abs=1e-6)
        with pytest.raises(ValueError, match="DEM"):
            tropo.stratified_tropo(delay, h[:10])
        with pytest.raises(ValueError, match=">= 50 pixels"):
            tropo.stratified_tropo(delay, h, np.zeros(h.shape, bool))
        with pytest.raises(ValueError, match="both x and y"):
            tropo.stratified_tropo(delay, h, x=np.arange(NX))
        with pytest.raises(ValueError, match="do not match"):
            tropo.stratified_tropo(delay, h, x=np.arange(5), y=np.arange(NY))


# --------------------------------------------------------------------------
# mode selection


def _dem_with_relief(relief_m: float, n: int = 64) -> np.ndarray:
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    # linear ramp: p95 - p5 of a uniform distribution is 0.9 of the range
    return relief_m / 0.9 * xx / (n - 1)


class TestChooseMode:
    def test_relief_is_p95_minus_p5(self):
        d = _dem_with_relief(900.0)
        assert tropo.dem_relief(d) == pytest.approx(900.0, rel=0.02)
        d[0, 0] = 1e5  # one tower does not count
        assert tropo.dem_relief(d) == pytest.approx(900.0, rel=0.02)
        d[:, :10] = np.nan
        assert np.isfinite(tropo.dem_relief(d))
        with pytest.raises(ValueError, match="two finite"):
            tropo.dem_relief(np.full((3, 3), np.nan))

    @pytest.mark.parametrize(
        ("relief", "expected"),
        [(130.0, "off"), (270.0, "off"), (410.0, "off"), (900.0, "off"),
         (1500.0, "stratified"), (2500.0, "stratified")],
    )
    def test_auto_bands(self, relief, expected):
        o = TropoOptions(mode="auto")
        assert tropo.choose_tropo_mode(o, _dem_with_relief(relief)) == expected

    def test_auto_needs_dem_and_bands_are_configurable(self):
        with pytest.raises(ValueError, match="needs the DEM"):
            tropo.choose_tropo_mode(TropoOptions(mode="auto"))
        o = TropoOptions(mode="auto", relief_off_meters=300, relief_stratified_meters=1000)
        assert tropo.choose_tropo_mode(o, _dem_with_relief(1200.0)) == "stratified"
        with pytest.raises(ValueError, match="must not exceed"):
            TropoOptions(mode="auto", relief_off_meters=2000, relief_stratified_meters=1000)

    def test_legacy_defers_to_apply_tropo_correction(self):
        o = TropoOptions()  # legacy
        assert tropo.choose_tropo_mode(o, apply_tropo_correction=True) == "full"
        assert tropo.choose_tropo_mode(o, apply_tropo_correction=False) == "off"
        for fixed in ("off", "stratified", "full"):
            assert tropo.choose_tropo_mode(TropoOptions(mode=fixed)) == fixed

    def test_frame_table_pins_the_mode(self):
        """T34.3: the per-frame entry wins over the relief rule."""
        table = load_frame_table()
        base = AlgorithmParameters()
        base.calibration_options.tropo.mode = "auto"
        la = table.apply(base, 16940).calibration_options
        assert la.tropo.mode == "stratified"
        assert tropo.choose_tropo_mode(la.tropo, _dem_with_relief(100.0)) == "stratified"
        houston = table.apply(base, 8882).calibration_options
        assert houston.tropo.mode == "off"
        assert tropo.choose_tropo_mode(houston.tropo, _dem_with_relief(3000.0)) == "off"
        # a frame without an entry keeps auto and follows the relief
        other = table.apply(base, 99999).calibration_options
        assert other.tropo.mode == "auto"
        assert tropo.choose_tropo_mode(other.tropo, _dem_with_relief(2000.0)) == "stratified"


class TestApplyTropo:
    def test_modes_produce_the_cal_tropo_component(self, dem):
        h = dem.values.astype(float)
        x, y = _grid_xy()
        rng = np.random.default_rng(4)
        corr = (-0.007 * h / 1e3 + 0.003 * rng.standard_normal(h.shape)).astype(
            np.float32
        )
        assert tropo.apply_tropo(corr, "off") is None
        assert tropo.apply_tropo(None, "off") is None
        full = tropo.apply_tropo(corr, "full")
        np.testing.assert_array_equal(full, corr)
        strat = tropo.apply_tropo(corr, "stratified", dem=h, x=x, y=y)
        assert strat is not None
        assert strat.dtype == np.float32
        ok = np.isfinite(h)
        # smoother than the input and carrying the height gradient
        assert np.nanstd(strat[ok] - (-0.007 * h / 1e3)[ok]) < 0.0003
        with pytest.raises(ValueError, match="needs the tropo correction"):
            tropo.apply_tropo(None, "full")
        with pytest.raises(ValueError, match="needs the DEM"):
            tropo.apply_tropo(corr, "stratified")
        with pytest.raises(ValueError, match="unknown"):
            tropo.apply_tropo(corr, "bogus")  # type: ignore[arg-type]

    def test_component_flows_into_calibrate_pair(self, dem):
        """Integration: stratified cal_tropo is restored into the calibration (D6)."""
        h = np.nan_to_num(dem.values.astype(float), nan=50.0)
        x, y = _grid_xy()
        xx, yy = np.meshgrid(x, y)
        gnss = 0.002 * (xx - X0) / (NX * PIXEL_M)
        rng = np.random.default_rng(5)
        tropo_full = -0.007 * h / 1e3 + 0.002 * np.sin(xx / 800.0)
        disp = gnss + 0.004 + tropo_full + 0.0005 * rng.standard_normal(h.shape)
        valid = np.ones(h.shape, bool)
        comp = tropo.apply_tropo(tropo_full, "stratified", dem=h, x=x, y=y)
        o = CalibrationOptions(
            surface={"method": "loclin", "cutoff_wavelength_meters": 20_000.0,
                     "two_pass": False},
        )
        res = calibrate_pair(disp, gnss, valid, (128, 128), o, PIXEL_M, 0.02773, tropo=comp)
        np.testing.assert_array_equal(res.cal_tropo, comp)
        assert res.components_applied["cal_tropo"]
        res.assert_closed()
        # the stratified part is gone from the product; the 5 km lateral part (2 pi x 800 m)
        # (no skill in the model, below the surface cutoff) stays in it
        calibrated = disp - res.calibration
        resid = calibrated - gnss
        assert np.corrcoef(resid.ravel(), (h / 1e3).ravel())[0, 1] ** 2 < 0.05
        assert np.corrcoef(resid.ravel(), np.sin(xx / 800.0).ravel())[0, 1] > 0.8


# --------------------------------------------------------------------------
# T34.4: independent recomputation on the golden pair (needs data + cal_disp)

GOLDEN = os.environ.get("CAL_DISP_GOLDEN_DIR")


@pytest.mark.skipif(not GOLDEN, reason="CAL_DISP_GOLDEN_DIR not set")
def test_golden_dem_surface_delay_matches_an_independent_interpolation():
    """Reproduces trade studies/tropo_check/independent_check.py: median within
    0.5 mm of a height-linear + bilinear lon/lat interpolation on the native
    cube at 3000 pixels (the trade study found median -0.02 mm, max 0.5 mm)."""
    cd = pytest.importorskip("cal_disp.product._tropo")
    golden = Path(GOLDEN)
    tropo_files = sorted((golden / "input_data" / "tropo").glob("*.nc"))
    dem_files = sorted((golden / "input_data" / "static_input").glob("*_dem.tif"))
    if not tropo_files or not dem_files:
        pytest.skip("golden dir has no tropo/DEM inputs")
    dem = rioxarray.open_rasterio(dem_files[0]).squeeze("band", drop=True)
    dem = dem.coarsen(x=6, y=6, boundary="trim").mean()  # 180 m, keeps the test short
    dem = dem.where(dem > -1000)  # the STATIC DEM uses a large negative nodata
    w, s, e, n = dem.rio.transform_bounds("EPSG:4326")
    prod = cd.TropoProduct.from_path(tropo_files[0])
    cube = prod.get_total_delay(
        bounds=(w, s, e, n), max_height=float(dem.max()) + 3e3, bounds_buffer=0.2
    )
    ours = tropo.interpolate_to_dem_surface(cube, dem).values.astype(float)

    h = dem.values.astype(float)
    rng = np.random.default_rng(0)
    iy = rng.integers(0, h.shape[0], 3000)
    ix = rng.integers(0, h.shape[1], 3000)
    ok = np.isfinite(h[iy, ix])
    iy, ix = iy[ok], ix[ok]
    lon, lat = Transformer.from_crs(dem.rio.crs, "EPSG:4326", always_xy=True).transform(
        dem.x.values[ix], dem.y.values[iy]
    )
    c = cube.transpose("height", "latitude", "longitude")
    latv, lonv, data = c.latitude.values, c.longitude.values, c.values
    if latv[0] > latv[-1]:
        latv, data = latv[::-1], data[:, ::-1, :]
    rgi = RegularGridInterpolator(
        (c.height.values, latv, lonv), data, bounds_error=False, fill_value=np.nan
    )
    independent = rgi(np.c_[h[iy, ix], lat, lon])
    d_mm = (ours[iy, ix] - independent) * 1e3
    assert np.isfinite(d_mm).mean() > 0.95
    assert abs(np.nanmedian(d_mm)) < 0.5
    assert np.nanpercentile(np.abs(d_mm), 99) < 2.0

