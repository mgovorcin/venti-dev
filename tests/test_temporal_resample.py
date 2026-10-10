# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Temporal resampling (PRD D15; plan T54)."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from venti.temporal import resample_dataarray, resample_timeseries

T0 = np.datetime64("2022-01-01")
NY, NX = 12, 15


def series(step_days=12, n=61, noise=0.003, seed=0, seasonal=0.0):
    """Linear + optional seasonal signal per pixel, white noise, 2 years at 12 days."""
    rng = np.random.default_rng(seed)
    dates = T0 + np.arange(n) * np.timedelta64(step_days, "D")
    t_years = np.arange(n) * step_days / 365.25
    _yy, xx = np.mgrid[0:NY, 0:NX].astype(float)
    rate = -0.02 + 0.01 * xx / NX  # m/yr
    truth = (
        rate[None] * t_years[:, None, None]
        + seasonal * np.sin(2 * np.pi * t_years)[:, None, None]
    )
    ts = truth + noise * rng.standard_normal((n, NY, NX))
    return ts, dates, truth, t_years


def truth_at(t_years_target, rate, seasonal):
    return (
        rate[None] * t_years_target[:, None, None]
        + seasonal * np.sin(2 * np.pi * t_years_target)[:, None, None]
    )


class TestLinear:
    @pytest.mark.parametrize("method", ["huber_linear", "linear", "theil_sen"])
    def test_denoises_a_linear_series_and_sigma_is_calibrated(self, method):
        ts, dates, _truth, _t_years = series(noise=0.004)
        targets = T0 + np.array([100, 203, 365, 500], dtype="timedelta64[D]")
        values, sigma = resample_timeseries(ts, dates, targets, 96.0, method)
        assert values.shape == sigma.shape == (4, NY, NX)
        rate = -0.02 + 0.01 * np.mgrid[0:NY, 0:NX][1] / NX
        exp = truth_at(np.array([100, 203, 365, 500]) / 365.25, rate, 0.0)
        err = values - exp
        # eight 12-day epochs per 96-day window: noise / sqrt(8) ~ 1.4 mm
        assert np.sqrt(np.nanmean(err**2)) < 0.0025
        assert abs(np.nanmean(err)) < 0.0005
        if method != "theil_sen":
            z = err / sigma
            assert 0.75 < np.nanstd(z) < 1.3, np.nanstd(z)

    def test_targets_on_acquisition_dates_match_the_truth_better_than_the_data(self):
        ts, dates, truth, _ = series(noise=0.005)
        values, _ = resample_timeseries(ts, dates, dates[10:14], 120.0)
        raw_err = np.sqrt(np.mean((ts[10:14] - truth[10:14]) ** 2))
        fit_err = np.sqrt(np.mean((values - truth[10:14]) ** 2))
        assert fit_err < 0.5 * raw_err


class TestWindow:
    def test_seasonal_bias_grows_with_the_window(self):
        """A short window follows the seasonal curve; a 1-year window flattens it."""
        ts, dates, _truth, _ = series(noise=0.001, seasonal=0.01)
        rate = -0.02 + 0.01 * np.mgrid[0:NY, 0:NX][1] / NX
        target_days = np.array([273, 365, 456])  # near the seasonal extremes
        targets = T0 + target_days.astype("timedelta64[D]")
        exp = truth_at(target_days / 365.25, rate, 0.01)
        errs = {}
        for window in (60.0, 120.0, 365.0):
            values, _ = resample_timeseries(ts, dates, targets, window, "linear")
            errs[window] = np.sqrt(np.nanmean((values - exp) ** 2))
        assert errs[60.0] < errs[120.0] < errs[365.0]
        assert errs[60.0] < 0.002  # seasonal amplitude is 10 mm
        assert errs[365.0] > 0.004

    def test_noise_falls_with_the_window_without_seasonal_signal(self):
        ts, dates, _truth, _ = series(noise=0.005, n=121, step_days=6)
        rate = -0.02 + 0.01 * np.mgrid[0:NY, 0:NX][1] / NX
        targets = T0 + np.array([200, 300, 400], dtype="timedelta64[D]")
        exp = truth_at(np.array([200, 300, 400]) / 365.25, rate, 0.0)
        e36, _ = resample_timeseries(ts, dates, targets, 36.0, "linear")
        e144, s144 = resample_timeseries(ts, dates, targets, 144.0, "linear")
        assert np.sqrt(np.mean((e144 - exp) ** 2)) < 0.6 * np.sqrt(
            np.mean((e36 - exp) ** 2)
        )
        # 24 epochs in the long window: sigma about 5 / sqrt(24) ~ 1 mm
        assert 0.0006 < np.nanmedian(s144) < 0.0016


class TestEdges:
    def test_series_ends_and_extrapolation(self):
        ts, dates, _truth, _ = series(noise=0.002)
        rate = -0.02 + 0.01 * np.mgrid[0:NY, 0:NX][1] / NX
        # inside the span but one-sided windows at both ends
        edge_days = np.array([0, 6, 714, 720])
        targets = T0 + edge_days.astype("timedelta64[D]")
        values, sigma = resample_timeseries(ts, dates, targets, 96.0)
        exp = truth_at(edge_days / 365.25, rate, 0.0)
        assert np.isfinite(values).all()
        assert np.sqrt(np.mean((values - exp) ** 2)) < 0.003
        # the one-sided fit is less certain than an interior one
        _, s_mid = resample_timeseries(ts, dates, T0 + np.timedelta64(360, "D"), 96.0)
        assert np.nanmedian(sigma) > 1.2 * np.nanmedian(s_mid)
        # outside the span: NaN unless extrapolation is allowed
        out = T0 + np.array([-10, 740], dtype="timedelta64[D]")
        v_out, _ = resample_timeseries(ts, dates, out, 96.0)
        assert np.isnan(v_out).all()
        v_ext, _ = resample_timeseries(ts, dates, out, 96.0, max_extrapolation_days=30)
        assert np.isfinite(v_ext).all()
        assert (
            np.sqrt(
                np.mean(
                    (v_ext - truth_at(np.array([-10, 740]) / 365.25, rate, 0.0)) ** 2
                )
            )
            < 0.004
        )

    def test_missing_epochs_per_pixel_and_min_points(self):
        ts, dates, _truth, _ = series(noise=0.002)
        ts[20:40, 0, 0] = np.nan  # a 240-day gap at one pixel
        ts[:, 1, 1] = np.nan  # a pixel with no data at all
        targets = T0 + np.array([360, 100], dtype="timedelta64[D]")
        values, sigma = resample_timeseries(ts, dates, targets, 96.0)
        assert np.isnan(values[0, 0, 0])
        assert np.isfinite(values[1, 0, 0])
        assert np.isnan(values[:, 1, 1]).all()
        assert np.isnan(sigma[:, 1, 1]).all()
        assert np.isfinite(values[:, 2:, 2:]).all()
        # a window with two epochs is refused by default
        ts2, dates2, *_ = series(step_days=40, n=10)
        v, _ = resample_timeseries(ts2, dates2, dates2[5], 90.0)  # 2-3 epochs
        v3, _ = resample_timeseries(ts2, dates2, dates2[5], 90.0, min_points=2)
        assert np.isnan(v).all() or np.isfinite(v3).all()

    def test_theil_sen_ignores_a_blundered_epoch(self):
        ts, dates, truth, _ = series(noise=0.001)
        ts[30] += 0.08  # one epoch 8 cm off (an unwrapping blunder)
        target = dates[30]
        v_ls, _ = resample_timeseries(ts, dates, target, 96.0, "linear")
        v_hub, _ = resample_timeseries(ts, dates, target, 96.0, "huber_linear")
        v_ts, _ = resample_timeseries(ts, dates, target, 96.0, "theil_sen")
        e = lambda v: np.sqrt(np.mean((v[0] - truth[30]) ** 2))  # noqa: E731
        assert e(v_ls) > 0.008  # 8 cm / 8 epochs
        assert e(v_hub) < 0.5 * e(v_ls)
        assert e(v_ts) < 0.002

    def test_input_checks_and_chunking(self):
        ts, dates, *_ = series()
        with pytest.raises(ValueError, match="epochs"):
            resample_timeseries(ts[:-1], dates, dates[0], 96.0)
        with pytest.raises(ValueError, match="window_days"):
            resample_timeseries(ts, dates, dates[0], 0.0)
        with pytest.raises(ValueError, match="unknown method"):
            resample_timeseries(ts, dates, dates[0], 96.0, "spline")  # type: ignore[arg-type]
        a, sa = resample_timeseries(ts, dates, dates[5:8], 96.0, chunk_pixels=7)
        b, sb = resample_timeseries(ts, dates, dates[5:8], 96.0)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(sa, sb)
        # float days and a 1-D series also work
        v, s = resample_timeseries(
            ts[:, 0, 0], np.arange(ts.shape[0]) * 12.0, [120.0], 96.0
        )
        assert v.shape == s.shape == (1,)


def test_dataarray_wrapper_keeps_coordinates():
    ts, dates, _truth, _ = series(noise=0.002)
    da = xr.DataArray(
        ts,
        dims=("time", "y", "x"),
        coords={"time": dates, "y": np.arange(NY) * 30.0, "x": np.arange(NX) * 30.0},
    )
    targets = T0 + np.array([120, 240], dtype="timedelta64[D]")
    out = resample_dataarray(da, targets, 96.0)
    assert set(out.data_vars) == {"value", "sigma"}
    assert out.value.dims == ("time", "y", "x")
    np.testing.assert_array_equal(out.time.values, targets)
    np.testing.assert_array_equal(out.y.values, da.y.values)
    v, s = resample_timeseries(ts, dates, targets, 96.0)
    np.testing.assert_allclose(out.value.values, v)
    np.testing.assert_allclose(out.sigma.values, s)
    assert out.attrs["window_days"] == 96.0
    with pytest.raises(ValueError, match="dimension"):
        resample_dataarray(da.rename({"time": "t"}), targets, 96.0)
    dask = pytest.importorskip("dask")  # noqa: F841 - research env only
    lazy = resample_dataarray(da.chunk({"y": 5, "x": 8}), targets, 96.0)
    np.testing.assert_allclose(lazy.value.values, v)
