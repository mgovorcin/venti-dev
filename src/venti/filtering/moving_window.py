"""Moving-window plane fitting for InSAR calibration.

Implements a moving-window least-squares surface fit to remove the
long-wavelength bias between InSAR displacement and a GNSS reference surface.
"""

from __future__ import annotations

import logging
import shutil
import tempfile

import numpy as np
from joblib import Parallel, delayed

from ..solver.plane_fitting import _fit_plane

logger = logging.getLogger(__name__)

# The plane is fit to every `_MAX_DECIMATE`-th pixel of a window for speed,
# but never to fewer than `_MIN_FIT_SAMPLES` pixels: with small windows or
# strong downsampling a fixed step left too few points and the fit silently
# returned a zero plane.
_MAX_DECIMATE = 50
_MIN_FIT_SAMPLES = 100


def _fit_decimation(shape: tuple[int, int]) -> int:
    """Largest step <= `_MAX_DECIMATE` keeping >= `_MIN_FIT_SAMPLES` pixels."""
    step = _MAX_DECIMATE
    while (
        step > 1 and (-(-shape[0] // step)) * (-(-shape[1] // step)) < _MIN_FIT_SAMPLES
    ):
        step //= 2
    return step


def _find_data_extent(array: np.ndarray, axis: int = 0) -> tuple[int, int]:
    """Return the start and stop indices of non-zero data along an axis.

    Parameters
    ----------
    array : np.ndarray
        Input array.
    axis : int, optional
        Axis to collapse before counting non-zeros, by default ``0``.

    Returns
    -------
    tuple of int
        ``(start, stop)`` indices.

    """
    count = np.count_nonzero(array, axis=axis)
    nonzero = np.where(count != 0)[0]
    if nonzero.size == 0:
        return 0, count.size
    return int(nonzero.min()), int(nonzero.max())


def _get_sliding_windows(
    length: int,
    win_size: int,
    win_overlap: int,
    first: int = 0,
    end: int = 0,
) -> list[slice]:
    """Build a list of overlapping 1-D window slices.

    Parameters
    ----------
    length : int
        Total length of the dimension.
    win_size : int
        Window size in pixels.
    win_overlap : int
        Overlap between adjacent windows in pixels.
    first : int, optional
        Start index, by default ``0``.
    end : int, optional
        End index; ``0`` means use `length`, by default ``0``.

    Returns
    -------
    list of slice
        Window slices covering ``[first, end)``.

    """
    if end == 0:
        end = length

    windows = []
    ix = 1
    stop = win_size
    while stop < end:
        start = first + (ix - 1) * win_size - win_overlap * (ix - 1)
        stop = first + ix * win_size - win_overlap * (ix - 1)
        ix += 1
        windows.append(slice(start, stop))

    if not windows:
        windows.append(slice(first, end))
    else:
        windows[-1] = slice(windows[-1].start, end)
    return windows


def _extend_window(
    slice_y: slice,
    slice_x: slice,
    extend_y: int,
    extend_x: int,
    length: int,
    width: int,
) -> tuple[tuple[slice, slice], tuple[slice, slice]]:
    """Expand a window by padding and return the extended slice and inner slice.

    Parameters
    ----------
    slice_y : slice
        Row slice of the original window.
    slice_x : slice
        Column slice of the original window.
    extend_y : int
        Pixels to extend in the row direction.
    extend_x : int
        Pixels to extend in the column direction.
    length : int
        Total number of rows.
    width : int
        Total number of columns.

    Returns
    -------
    extended_win : tuple of slice
        ``(row_slice, col_slice)`` for the extended window.
    padding : tuple of slice
        ``(row_slice, col_slice)`` indexing the original window inside the
        extended window.

    """
    y_start = max(slice_y.start - extend_y, 0)
    y_stop = min(slice_y.stop + extend_y, length)
    x_start = max(slice_x.start - extend_x, 0)
    x_stop = min(slice_x.stop + extend_x, width)

    extended_win = (slice(y_start, y_stop), slice(x_start, x_stop))
    nlength = y_stop - y_start
    nwidth = x_stop - x_start
    padding = (
        slice(slice_y.start - y_start, nlength - (y_stop - slice_y.stop)),
        slice(slice_x.start - x_start, nwidth - (x_stop - slice_x.stop)),
    )
    return extended_win, padding


def _get_coordinate_grid(
    win_y: list[int],
    win_x: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Build a normalised coordinate meshgrid for a window region.

    Coordinates are expressed as normalised values in ``[-1, 1]`` relative to
    the window extent rather than in absolute geographic or projected units.
    This keeps the polynomial design matrix well-conditioned regardless of
    whether the bounds are in lat/lon degrees or UTM metres.

    Parameters
    ----------
    win_y : list of int
        ``[y_start, y_stop]`` pixel row indices of the window.
    win_x : list of int
        ``[x_start, x_stop]`` pixel column indices of the window.

    Returns
    -------
    grid_x : np.ndarray
        Normalised x (column) coordinates, shape ``(ny, nx)``, range ``[-1, 1]``.
    grid_y : np.ndarray
        Normalised y (row) coordinates, shape ``(ny, nx)``, range ``[-1, 1]``.

    """
    nx = win_x[1] - win_x[0]
    ny = win_y[1] - win_y[0]
    xs = np.linspace(-1.0, 1.0, nx)
    ys = np.linspace(-1.0, 1.0, ny)
    return np.meshgrid(xs, ys)


def _process_window(
    win_index: tuple[slice, slice],
    insar_data: np.ndarray,
    gnss_los: np.ndarray,
    win_extend_y: int,
    win_extend_x: int,
    length: int,
    width: int,
    poly_order: float,
    gnss_los_std: np.ndarray | None = None,
) -> tuple[
    tuple[slice, slice], np.ndarray | None, np.ndarray | None, np.ndarray | None
]:
    """Fit a calibration plane within a single moving window.

    Designed to be called in parallel via ``joblib``.

    Parameters
    ----------
    win_index : tuple of slice
        ``(row_slice, col_slice)`` of the window to fill.
    insar_data : np.ndarray
        Full InSAR data array (gap-filled, masked).
    gnss_los : np.ndarray
        Full GNSS LOS reference array.
    win_extend_y : int
        Row extension for fitting context.
    win_extend_x : int
        Column extension for fitting context.
    length : int
        Total number of rows.
    width : int
        Total number of columns.
    poly_order : float
        Polynomial order for plane fitting.
    gnss_los_std : np.ndarray, optional
        Full GNSS LOS uncertainty array, same shape as `gnss_los`. When
        given, the window's slice of it is passed to `_fit_plane` as
        inverse-variance weights, so pixels near a well-constrained GNSS
        point dominate the local fit more than pixels far from any station
        (where the interpolated reference is least trustworthy). ``None``
        (default) fits with uniform weights.

    Returns
    -------
    win_index : tuple of slice
    plane : np.ndarray or None
        Plane over `win_index`, 0 where the data are invalid.
    plane_std : np.ndarray or None
    valid : np.ndarray or None
        Where the data are valid; only these pixels may carry weight when
        windows are blended.

    """
    win2, pad = _extend_window(
        win_index[0], win_index[1], win_extend_y, win_extend_x, length, width
    )
    win_y = [win2[0].start, win2[0].stop]
    win_x = [win2[1].start, win2[1].stop]
    win_lons, win_lats = _get_coordinate_grid(win_y=win_y, win_x=win_x)

    res = insar_data[win2] - gnss_los[win2]
    if np.isnan(res).sum() / res.size > 0.8:
        return win_index, None, None, None

    res_std = gnss_los_std[win2] if gnss_los_std is not None else None
    try:
        plane, plane_std = _fit_plane(
            res,
            win_lons,
            win_lats,
            order=poly_order,
            decimate=_fit_decimation(res.shape),
            data_std=res_std,
        )
        valid = np.isfinite(res[pad])
        filled_plane = np.where(valid, plane[pad], 0.0)
        filled_std = np.where(valid, plane_std[pad], 0.0)
    except Exception:
        logger.debug("Skipping window %s", win_index, exc_info=True)
        return win_index, None, None, None
    else:
        return win_index, filled_plane, filled_std, valid


def fit_windowed_plane(
    insar_data: np.ndarray,
    gnss_los: np.ndarray,
    win_xsize: int,
    win_ysize: int,
    win_overlap_x: int,
    win_overlap_y: int,
    win_extend_x: int,
    win_extend_y: int,
    gnss_los_std: np.ndarray | None = None,
    poly_order: float = 1.5,
    n_jobs: int = -1,
    smoothing_sigma: float | None = None,
    smoothing_method: str = "gaussian",
    sg_window_length: int = 51,
    sg_polyorder: int = 3,
    mask_residual_outliers: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a windowed polynomial calibration surface to InSAR data.

    Uses a moving-window least-squares approach to estimate the
    long-wavelength surface ``m`` such that ``insar_data ≈ gnss_los + m``.
    The correction is applied as ``calibrated = insar_data - m``.

    Parameters
    ----------
    insar_data : np.ndarray
        2-D InSAR displacement or velocity array.
    gnss_los : np.ndarray
        2-D GNSS LOS reference field (same shape as `insar_data`).
    win_xsize : int
        Window width in pixels.
    win_ysize : int
        Window height in pixels.
    win_overlap_x : int
        Window overlap in the x direction.
    win_overlap_y : int
        Window overlap in the y direction.
    win_extend_x : int
        Extension beyond the window for fitting context.
    win_extend_y : int
        Extension beyond the window for fitting context.
    gnss_los_std : np.ndarray, optional
        Uncertainty of the GNSS LOS field, same shape as `gnss_los`. When
        given, used as inverse-variance weights in each window's
        least-squares plane fit (see `_process_window`/`_fit_plane`), so
        pixels near a well-constrained GNSS point influence the local fit
        more than pixels far from any station. ``None`` (default) fits
        with uniform weights.
    poly_order : float, optional
        Polynomial order for plane fitting, by default ``1.5``.
    n_jobs : int, optional
        Number of parallel jobs for ``joblib.Parallel``,
        by default ``-1`` (all CPUs).
    smoothing_sigma : float, optional
        Standard deviation (pixels) passed to the post-assembly low-pass
        filter.  Ignored when ``smoothing_method="savitzky_golay"``.
        ``None`` disables smoothing.
    smoothing_method : str, optional
        Post-assembly low-pass filter to apply.  One of ``"gaussian"``
        (default, spatial-domain), ``"gaussian_fft"``, ``"hanning_fft"``,
        or ``"savitzky_golay"``.
    sg_window_length : int, optional
        Window length for the Savitzky-Golay filter in pixels (must be odd),
        by default ``51``.  Only used when ``smoothing_method="savitzky_golay"``.
    sg_polyorder : int, optional
        Polynomial order for the Savitzky-Golay filter, by default ``3``.
        Only used when ``smoothing_method="savitzky_golay"``.
    mask_residual_outliers : bool, optional
        Whether to exclude the most extreme 15% / 85% of
        ``insar_data - gnss_los`` residual pixels (per `_get_residual_mask`)
        before filling and fitting, by default ``True``. This is a fixed
        quantile trim, unconditional and independent of any caller-level
        outlier detection (e.g. `process_displacement_file`'s
        ``residual_outlier_mad_threshold``) — it always discards ~30% of
        pixels, real deformation included, and replaces them with
        InSAR-only smoothed values before the fit ever sees them. Set to
        ``False`` for a fit that uses `insar_data` as given, masking only
        genuinely invalid (NaN) pixels.

    Returns
    -------
    calibration_surface : np.ndarray
        Estimated calibration surface to subtract from `insar_data`.
    calibration_std : np.ndarray
        Uncertainty of the calibration surface.

    Examples
    --------
    ::

        surface, surface_std = fit_windowed_plane(
            insar_data=displacement_mm,
            gnss_los=gnss_los_mm,
            win_xsize=1000, win_ysize=1000,
            win_overlap_x=10, win_overlap_y=10,
            win_extend_x=1000, win_extend_y=1000,
        )
        calibrated = displacement_mm - surface

    """
    from ..spatial.gap_filling import _get_residual_mask, fill_gaps

    invalid_mask = np.ma.masked_invalid(insar_data).mask
    outlier_mask = (
        _get_residual_mask(insar_data, gnss_los)
        if mask_residual_outliers
        else invalid_mask
    )
    length, width = insar_data.shape

    logger.info(
        "Fitting windowed calibration plane: %d x %d px windows, poly_order=%s",
        win_ysize,
        win_xsize,
        poly_order,
    )

    insar_filled = fill_gaps(
        insar_data,
        outlier_mask,
        max_search_distance=max(insar_data.shape) // 4,
        smoothing_iterations=10,
    )
    insar_filled = np.where(invalid_mask, np.nan, insar_filled)

    y_start, y_stop = _find_data_extent(insar_data, axis=1)
    win_ys = _get_sliding_windows(length, win_ysize, win_overlap_y, y_start, y_stop)

    all_windows: list[tuple[slice, slice]] = []
    for win_y in win_ys:
        x_start, x_stop = _find_data_extent(insar_data[win_y, :], axis=0)
        win_xs = _get_sliding_windows(
            x_stop - x_start, win_xsize, win_overlap_x, x_start, x_stop
        )
        all_windows.extend((win_y, win_x) for win_x in win_xs)

    logger.info("Processing %d windows (n_jobs=%d)", len(all_windows), n_jobs)

    _tmpdir = tempfile.mkdtemp(prefix="venti_fit_")
    try:
        _insar_path = f"{_tmpdir}/insar.mmap"
        _gnss_path = f"{_tmpdir}/gnss.mmap"

        _insar_mm = np.memmap(
            _insar_path, dtype=insar_filled.dtype, mode="w+", shape=insar_filled.shape
        )
        _gnss_mm = np.memmap(
            _gnss_path, dtype=gnss_los.dtype, mode="w+", shape=gnss_los.shape
        )
        _insar_mm[:] = insar_filled
        _gnss_mm[:] = gnss_los
        del _insar_mm, _gnss_mm

        _insar_mm = np.memmap(
            _insar_path, dtype=insar_filled.dtype, mode="r", shape=insar_filled.shape
        )
        _gnss_mm = np.memmap(
            _gnss_path, dtype=gnss_los.dtype, mode="r", shape=gnss_los.shape
        )

        _gnss_std_mm = None
        if gnss_los_std is not None:
            _gnss_std_path = f"{_tmpdir}/gnss_std.mmap"
            _gnss_std_mm = np.memmap(
                _gnss_std_path,
                dtype=gnss_los_std.dtype,
                mode="w+",
                shape=gnss_los_std.shape,
            )
            _gnss_std_mm[:] = gnss_los_std
            del _gnss_std_mm
            _gnss_std_mm = np.memmap(
                _gnss_std_path,
                dtype=gnss_los_std.dtype,
                mode="r",
                shape=gnss_los_std.shape,
            )

        results = Parallel(n_jobs=n_jobs)(
            delayed(_process_window)(
                ix,
                _insar_mm,
                _gnss_mm,
                win_extend_y,
                win_extend_x,
                length,
                width,
                poly_order,
                gnss_los_std=_gnss_std_mm,
            )
            for ix in all_windows
        )
    finally:
        shutil.rmtree(_tmpdir, ignore_errors=True)

    cal_surface = np.zeros(insar_data.shape, dtype=np.float64)
    cal_std = np.zeros(insar_data.shape, dtype=np.float64)
    weight_sum = np.zeros(insar_data.shape, dtype=np.float64)
    weight_std = np.zeros(insar_data.shape, dtype=np.float64)

    for ix, plane_val, std_val, valid in results:
        if plane_val is None:
            continue
        ny, nx = plane_val.shape
        # 2-D Hann taper: higher weight toward window centre, tapers to zero at edges.
        # Invalid pixels (water, outside the swath) get no weight: counting their
        # 0s would pull the surface toward 0 near every mask edge (and the
        # smoothing spreads that inland), so the fit would no longer be linear
        # in `gnss_los` there.
        taper = np.outer(np.hanning(ny), np.hanning(nx)) * valid
        cal_surface[ix] += plane_val * taper
        weight_sum[ix] += taper
        if std_val is not None:
            cal_std[ix] += std_val * taper
            weight_std[ix] += taper

    nonzero = weight_sum > 0
    cal_surface[nonzero] /= weight_sum[nonzero]
    nonzero_std = weight_std > 0
    cal_std[nonzero_std] /= weight_std[nonzero_std]

    _apply_smoothing = (
        smoothing_sigma is not None or smoothing_method == "savitzky_golay"
    )
    if _apply_smoothing:
        from .low_pass_filters import (
            gaussian_fft,
            gaussian_spatial,
            hanning_fft,
            savitzky_golay,
        )

        valid = nonzero
        if smoothing_method == "gaussian_fft":
            cal_surface = gaussian_fft(cal_surface, valid, smoothing_sigma)  # type: ignore[arg-type]
            if nonzero_std.any():
                cal_std = gaussian_fft(cal_std, nonzero_std, smoothing_sigma)  # type: ignore[arg-type]
        elif smoothing_method == "hanning_fft":
            cal_surface = hanning_fft(cal_surface, valid, smoothing_sigma)  # type: ignore[arg-type]
            if nonzero_std.any():
                cal_std = hanning_fft(cal_std, nonzero_std, smoothing_sigma)  # type: ignore[arg-type]
        elif smoothing_method == "savitzky_golay":
            cal_surface = savitzky_golay(
                cal_surface, valid, sg_window_length, sg_polyorder
            )
            if nonzero_std.any():
                cal_std = savitzky_golay(
                    cal_std, nonzero_std, sg_window_length, sg_polyorder
                )
        else:
            # default: spatial Gaussian
            cal_surface = gaussian_spatial(cal_surface, valid, smoothing_sigma)  # type: ignore[arg-type]
            if nonzero_std.any():
                cal_std = gaussian_spatial(cal_std, nonzero_std, smoothing_sigma)  # type: ignore[arg-type]

    return cal_surface.astype(np.float32), cal_std.astype(np.float32)
