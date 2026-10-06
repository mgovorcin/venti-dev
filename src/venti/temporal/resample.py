# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Moving-window denoising and temporal interpolation (PRD D15; plan T54).

Ascending and descending stacks are acquired on different dates. To
decompose them per pair (`venti.decomposition`) both are resampled onto
common epochs: for each target epoch ``t0`` and each pixel, a robust linear
fit to the calibrated displacement within ``|t - t0| <= window_days / 2`` is
evaluated at ``t0``. The fit denoises (the window averages the DISP noise)
and interpolates (the line crosses the gaps between acquisitions); its
uncertainty is the standard error of the fitted line at ``t0`` from the
window residuals, ``sigma_r * sqrt(1 / n_eff + (t0 - t_bar)^2 / S_tt)``.

Methods:
- ``"huber_linear"`` (default): weighted least squares with Huber
  re-weighting (two passes) of residuals scaled by 1.4826 MAD; vectorised
  over pixels through the 2x2 normal equations;
- ``"theil_sen"``: median of pairwise slopes, median intercept; immune to a
  single blundered epoch but O(n^2) in the window length;
- ``"linear"``: plain least squares.

The window should be short against the seasonal period (a 1-year window
fits a line through a full seasonal cycle and biases the value at ``t0`` by
the curvature); 60-120 days for Sentinel-1 6/12-day sampling gives 5-20
epochs per window. The series ends are one-sided windows: the value is
still a fit at ``t0`` as long as ``t0`` lies inside the data span (plus
`max_extrapolation_days`), with a larger sigma because ``t0`` is far from
``t_bar``. Spatial chunking keeps the memory at ``chunk_pixels x window``.

"""

from __future__ import annotations

import logging
import warnings
from datetime import datetime
from typing import Any, Literal

import numpy as np

logger = logging.getLogger(__name__)

__all__ = ["ResampleMethod", "resample_dataarray", "resample_timeseries"]

ResampleMethod = Literal["huber_linear", "theil_sen", "linear"]
MAD_SCALE = 1.4826
HUBER_C = 1.345
EPS = 1e-12


def _to_days(dates: Any) -> np.ndarray:
    """Convert datetimes / numpy datetime64 / numbers to float days."""
    arr = np.asarray(dates)
    if np.issubdtype(arr.dtype, np.datetime64):
        return (arr.astype("datetime64[ns]").astype(np.int64) / 86_400e9).astype(
            np.float64
        )
    if arr.dtype == object and arr.size and isinstance(arr.flat[0], datetime):
        return np.array(
            [np.datetime64(d, "ns").astype(np.int64) / 86_400e9 for d in arr.ravel()]
        ).reshape(arr.shape)
    return arr.astype(np.float64)


def _weighted_line(
    t: np.ndarray, y: np.ndarray, w: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Weighted LS line per pixel. t: (n,), y/w: (n, p). Returns a, b, t_bar, S_tt."""
    s0 = w.sum(0)
    s0_safe = np.where(s0 > 0, s0, 1.0)
    t_bar = (w * t[:, None]).sum(0) / s0_safe
    dt = t[:, None] - t_bar
    s_tt = (w * dt * dt).sum(0)
    y_bar = (w * y).sum(0) / s0_safe
    s_ty = (w * dt * (y - y_bar)).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        b = np.where(s_tt > EPS, s_ty / np.where(s_tt > EPS, s_tt, 1.0), np.nan)
    a = y_bar  # value at t_bar
    return a, b, t_bar, s_tt


def _theil_sen(
    t: np.ndarray, y: np.ndarray, ok: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Theil-Sen slope and intercept per pixel (y: (n, p), ok: (n, p) finite mask)."""
    n = t.size
    i, j = np.triu_indices(n, k=1)
    dt = (t[j] - t[i])[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        slopes = (y[j] - y[i]) / dt
    pair_ok = ok[i] & ok[j] & (np.abs(dt) > EPS)
    slopes = np.where(pair_ok, slopes, np.nan)
    b = np.nanmedian(slopes, axis=0)
    inter = np.where(ok, y - b[None, :] * t[:, None], np.nan)
    a0 = np.nanmedian(inter, axis=0)
    return a0, b


def _fit_window(
    t: np.ndarray,
    y: np.ndarray,
    t0: float,
    method: ResampleMethod,
    min_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit each pixel's window (y: (n, p)) and evaluate at t0 -> value, sigma, n_eff."""
    ok = np.isfinite(y)
    w = ok.astype(np.float64)
    n_ok = ok.sum(0)
    y0 = np.where(ok, y, 0.0)
    with warnings.catch_warnings():
        # pixels with no finite epoch in the window give all-NaN medians
        warnings.simplefilter("ignore", RuntimeWarning)
        return _fit_window_impl(t, y, y0, ok, w, n_ok, t0, method, min_points)


def _fit_window_impl(
    t: np.ndarray,
    y: np.ndarray,
    y0: np.ndarray,
    ok: np.ndarray,
    w: np.ndarray,
    n_ok: np.ndarray,
    t0: float,
    method: ResampleMethod,
    min_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if method == "theil_sen":
        a0, b = _theil_sen(t, np.where(ok, y, np.nan), ok)
        value = a0 + b * t0
        resid = np.where(ok, y0 - (a0 + b * t[:, None]), 0.0)
        t_bar = (w * t[:, None]).sum(0) / np.maximum(n_ok, 1)
        s_tt = (w * (t[:, None] - t_bar) ** 2).sum(0)
        n_eff = n_ok.astype(np.float64)
    else:
        a, b, t_bar, s_tt = _weighted_line(t, y0, w)
        if method == "huber_linear":
            for _ in range(2):
                resid = np.where(ok, y0 - (a + b * (t[:, None] - t_bar)), 0.0)
                med = np.nanmedian(np.where(ok, resid, np.nan), axis=0)
                mad = MAD_SCALE * np.nanmedian(
                    np.where(ok, np.abs(resid - med), np.nan), axis=0
                )
                scale = np.where(mad > EPS, mad, np.nan)
                with np.errstate(invalid="ignore", divide="ignore"):
                    z = np.abs(resid - med) / scale
                    hub = np.where(z > HUBER_C, HUBER_C / z, 1.0)
                hub = np.where(np.isfinite(hub), hub, 1.0)
                w = ok * hub
                a, b, t_bar, s_tt = _weighted_line(t, y0, w)
        value = a + b * (t0 - t_bar)
        resid = np.where(ok, y0 - (a + b * (t[:, None] - t_bar)), 0.0)
        s0 = w.sum(0)
        s2 = (w * w).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            n_eff = np.where(s2 > 0, s0 * s0 / s2, 0.0)
    # residual scale from every finite epoch (not the robust weights, which
    # would shrink the apparent scatter) with the 2-parameter dof, and the
    # standard error of the line at t0; a blundered epoch inflates sigma,
    # which is the conservative side
    with np.errstate(invalid="ignore", divide="ignore"):
        dof = np.maximum(n_ok - 2.0, 1.0)
        sigma_r = np.sqrt((resid**2).sum(0) / dof)
        lever = (t0 - t_bar) ** 2 / np.where(s_tt > EPS, s_tt, np.nan)
        sigma = sigma_r * np.sqrt(1.0 / np.maximum(n_eff, 1.0) + lever)
    bad = (n_ok < min_points) | ~np.isfinite(value) | (s_tt <= EPS)
    value = np.where(bad, np.nan, value)
    sigma = np.where(bad, np.nan, sigma)
    return value, sigma, np.where(bad, 0.0, n_eff)


def resample_timeseries(
    ts: np.ndarray,
    dates: Any,
    target_dates: Any,
    window_days: float,
    method: ResampleMethod = "huber_linear",
    *,
    min_points: int = 3,
    max_extrapolation_days: float = 0.0,
    chunk_pixels: int = 200_000,
) -> tuple[np.ndarray, np.ndarray]:
    """Resample a time series onto `target_dates` with moving-window line fits.

    Parameters
    ----------
    ts : np.ndarray
        Displacement with time first, ``(nt, ...)``; NaN where missing.
    dates : array-like
        Acquisition epochs (datetime64, datetime or float days), length nt.
    target_dates : array-like
        Epochs to evaluate at, length m.
    window_days : float
        Full window width; epochs with ``|t - t0| <= window_days / 2`` are fitted.
    method : {"huber_linear", "theil_sen", "linear"}
        Robust fit (see module docstring).
    min_points : int
        Fewer finite epochs in the window give NaN.
    max_extrapolation_days : float
        Targets farther than this outside the series span are NaN.
    chunk_pixels : int
        Pixels processed per block (memory ``chunk x window``).

    Returns
    -------
    values, sigma : np.ndarray
        ``(m, ...)`` arrays: the denoised value at each target and its
        standard error.

    """
    y_all = np.asarray(ts, dtype=np.float64)
    t = _to_days(dates)
    t_out = np.atleast_1d(_to_days(target_dates))
    if y_all.ndim < 1 or y_all.shape[0] != t.size:
        msg = f"ts has {y_all.shape[0] if y_all.ndim else 0} epochs, dates has {t.size}"
        raise ValueError(msg)
    if window_days <= 0 or min_points < 2 or chunk_pixels < 1:
        msg = "window_days must be > 0, min_points >= 2 and chunk_pixels >= 1"
        raise ValueError(msg)
    if method not in ("huber_linear", "theil_sen", "linear"):
        msg = f"unknown method {method!r}"
        raise ValueError(msg)
    space = y_all.shape[1:]
    p = int(np.prod(space)) if space else 1
    y2 = y_all.reshape(t.size, p)
    order = np.argsort(t, kind="stable")
    t, y2 = t[order], y2[order]
    t_lo, t_hi = t[0] - max_extrapolation_days, t[-1] + max_extrapolation_days
    half = window_days / 2.0

    values = np.full((t_out.size, p), np.nan)
    sigmas = np.full((t_out.size, p), np.nan)
    for k, t0 in enumerate(t_out):
        if not (t_lo <= t0 <= t_hi):
            continue
        sel = np.abs(t - t0) <= half
        if sel.sum() < min_points:
            continue
        tw, yw = t[sel], y2[sel]
        for start in range(0, p, chunk_pixels):
            block = slice(start, min(start + chunk_pixels, p))
            v, s, _ = _fit_window(tw, yw[:, block], float(t0), method, min_points)
            values[k, block] = v
            sigmas[k, block] = s
    n_valid = int(np.isfinite(values).sum())
    logger.debug(
        "resample_timeseries(%s, %g d): %d targets, %.1f%% valid",
        method,
        window_days,
        t_out.size,
        100.0 * n_valid / max(values.size, 1),
    )
    return values.reshape((t_out.size, *space)), sigmas.reshape((t_out.size, *space))


def resample_dataarray(
    da: Any,
    target_dates: Any,
    window_days: float,
    method: ResampleMethod = "huber_linear",
    time_dim: str = "time",
    **kwargs: Any,
) -> Any:
    """`resample_timeseries` on an ``xarray.DataArray``; returns a Dataset.

    The result has ``value`` and ``sigma`` on ``(time, <other dims>)`` with
    the target dates as the time coordinate. A dask-backed input is
    processed chunk by chunk over space through ``xarray.apply_ufunc``
    (``dask="parallelized"``); the time dimension must be a single chunk.
    """
    import xarray as xr

    if time_dim not in da.dims:
        msg = f"{time_dim!r} is not a dimension of the array {da.dims}"
        raise ValueError(msg)
    da = da.transpose(time_dim, ...)
    targets = np.atleast_1d(np.asarray(target_dates))
    dates = da[time_dim].values

    def _core(block: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # apply_ufunc moves core dims last: block is (..., nt)
        moved = np.moveaxis(block, -1, 0)
        v, s = resample_timeseries(moved, dates, targets, window_days, method, **kwargs)
        return np.moveaxis(v, 0, -1), np.moveaxis(s, 0, -1)

    is_dask = getattr(da.data, "chunks", None) is not None
    value, sigma = xr.apply_ufunc(
        _core,
        da,
        input_core_dims=[[time_dim]],
        output_core_dims=[["target_time"], ["target_time"]],
        dask="parallelized" if is_dask else "forbidden",
        output_dtypes=[np.float64, np.float64],
        dask_gufunc_kwargs=(
            {"output_sizes": {"target_time": targets.size}} if is_dask else None
        ),
    )
    out = xr.Dataset({"value": value, "sigma": sigma})
    out = out.assign_coords(target_time=targets).rename({"target_time": time_dim})
    out = out.transpose(time_dim, ...)
    out.attrs.update(
        {"window_days": window_days, "method": method, "source_time_dim": time_dim}
    )
    return out
