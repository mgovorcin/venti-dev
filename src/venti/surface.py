"""Estimate a GNSS calibration surface from arrays.

The numerical core of Venti's calibration, independent of files and
configuration objects, so it can be used in any workflow. The file-based
wrapper is `venti.workflow.calibration.process_displacement_file`.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Sequence
from contextlib import nullcontext
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import threading

    from .workflow.config import CalibrationOptions

logger = logging.getLogger(__name__)

SENTINEL1_WAVELENGTH_M = 0.05546


@dataclass(frozen=True)
class CalibrationSurface:
    """Result of `estimate_calibration_surface`.

    Attributes
    ----------
    surface : np.ndarray
        Calibration surface in meters, same shape and precision as the input
        displacement.
        ``disp - surface`` is the calibrated displacement.
    n_auto_masked_pixels : int or None
        Pixels filled by the MAD outlier test (downsampled grid), or ``None``
        if the test did not run.
    n_region_masked_pixels : int or None
        Pixels filled by the coherent-region test (downsampled grid), or
        ``None`` if the test did not run.

    """

    surface: np.ndarray
    n_auto_masked_pixels: int | None = None
    n_region_masked_pixels: int | None = None


def estimate_calibration_surface(
    disp: np.ndarray,
    gnss_los: np.ndarray,
    mask: np.ndarray,
    ref_point: tuple[int, int],
    window_size: int | tuple[int, int],
    *,
    corrections: Sequence[np.ndarray] = (),
    gnss_los_std: np.ndarray | None = None,
    event_mask: np.ndarray | None = None,
    options: CalibrationOptions | None = None,
    wavelength_m: float = SENTINEL1_WAVELENGTH_M,
    downsample_factor: int = 1,
    downsample_method: str = "mean",
    downsample_weights: np.ndarray | None = None,
    n_jobs: int = -1,
    fit_lock: threading.Lock | None = None,
    fit_surface: Callable[..., np.ndarray] | None = None,
) -> CalibrationSurface:
    """Estimate the surface that ties an InSAR displacement field to GNSS.

    Steps: remove `corrections`, zero at `ref_point`, mask, correct
    unwrapping errors, fill event regions, downsample, optionally detect
    unmasked deformation, fit windowed planes to ``disp - gnss_los``,
    upsample, and add back the corrections and reference offset.
    Corrections are removed before the fit because the fit is close to
    linear, ``fit(disp, gnss) = fit(disp - c, gnss) + fit(c, 0)``: leaving a
    smooth correction in would let the fit absorb part of it.

    Parameters
    ----------
    disp : np.ndarray
        LOS displacement in meters, 2-D.
    gnss_los : np.ndarray
        GNSS LOS displacement for the same interval, in meters, same shape.
    mask : np.ndarray
        Boolean valid-pixel mask, same shape.
    ref_point : tuple of int
        Reference pixel ``(row, col)``; InSAR has no absolute datum, so the
        displacement is zeroed there and the offset is added back.
    window_size : int or tuple of int
        Fit window size in pixels, ``n`` or ``(nx, ny)``, before downsampling.
    corrections : sequence of np.ndarray, optional
        Signals GNSS does not contain (e.g. tropospheric delay, solid Earth
        tide), in meters, each to be subtracted from `disp`.
    gnss_los_std : np.ndarray, optional
        GNSS LOS uncertainty in meters, used as inverse-variance fit weights.
    event_mask : np.ndarray, optional
        True = valid, False = deformation to exclude; excluded pixels are
        filled from their neighbours before the fit.
    options : CalibrationOptions, optional
        Algorithm options; defaults to ``CalibrationOptions()``.
    wavelength_m : float, optional
        Radar wavelength for unwrapping-error correction, by default
        Sentinel-1 C-band.
    downsample_factor : int, optional
        Fit on a grid downsampled by this factor, by default 1 (none).
    downsample_method : str, optional
        ``'mean'`` or ``'median'``, by default ``'mean'``.
    downsample_weights : np.ndarray, optional
        Per-pixel weights for downsampling (e.g. temporal coherence).
    n_jobs : int, optional
        Parallel workers for the windowed fit, by default -1 (all CPUs).
    fit_lock : threading.Lock, optional
        Held during the windowed fit only, to serialise fits across threads.
    fit_surface : callable, optional
        Replacement for `SpatialProcessor.fit_windowed_surface` (same
        keyword arguments), e.g. to test or swap the fitting method.

    Returns
    -------
    CalibrationSurface
        The surface (meters) and outlier-detection counts.

    Raises
    ------
    ValueError
        If an input array or correction does not match `disp`'s shape, or
        the displacement is NaN at `ref_point`.

    Examples
    --------
    >>> import numpy as np
    >>> from venti.surface import estimate_calibration_surface
    >>> from venti.workflow.config import CalibrationOptions
    >>> yy, xx = np.mgrid[:60, :80]
    >>> gnss_los = 0.002 * xx / 80
    >>> disp = gnss_los + 0.004 * yy / 60  # plus an orbital-like ramp
    >>> result = estimate_calibration_surface(
    ...     disp, gnss_los, np.ones(disp.shape, bool), (30, 40), 20,
    ...     options=CalibrationOptions(unwrap_error_correction=False), n_jobs=1,
    ... )
    >>> calibrated = disp - result.surface
    >>> result.surface.shape
    (60, 80)

    """
    from scipy.ndimage import binary_dilation

    from .spatial.interpolation import (
        detect_coherent_residual_regions,
        fill_masked_region,
    )
    from .spatial.resample import downsample_array, upsample_array
    from .unwrap import correct_region_offset
    from .workflow.config import CalibrationOptions

    opts = options if options is not None else CalibrationOptions()
    for name, arr in [
        ("gnss_los", gnss_los),
        ("mask", mask),
        ("gnss_los_std", gnss_los_std),
        ("event_mask", event_mask),
        ("downsample_weights", downsample_weights),
        *[(f"corrections[{i}]", c) for i, c in enumerate(corrections)],
    ]:
        if arr is not None and np.shape(arr) != disp.shape:
            msg = f"{name} shape {np.shape(arr)} != disp shape {disp.shape}"
            raise ValueError(msg)
    if isinstance(window_size, int | np.integer):
        win_x = win_y = int(window_size)
    else:
        win_x, win_y = window_size

    # Keep the input precision (float32 for OPERA products): at full
    # resolution each extra float64 copy of a frame costs ~0.6 GB.
    disp = np.array(disp, dtype=np.promote_types(disp.dtype, np.float32))
    total_correction = np.zeros(disp.shape, dtype=disp.dtype)
    for correction in corrections:
        total_correction += correction
    disp -= total_correction

    refy, refx = ref_point
    disp0 = float(disp[refy, refx])
    # A NaN here would turn the whole epoch into NaN. Alternative kept for
    # later: shift by the median of valid pixels, which the shift-equivariant
    # fit makes equivalent and which no single pixel can break.
    if not np.isfinite(disp0):
        msg = (
            f"Reference pixel {ref_point} has no valid displacement "
            f"(value {disp0}); choose a reference_point valid in every epoch"
        )
        raise ValueError(msg)
    disp -= disp0

    disp = np.where(mask & ~np.isnan(disp), disp, np.nan)

    if event_mask is not None:
        event_mask = np.asarray(event_mask, dtype=bool)
        if opts.event_mask_buffer_pixels > 0:
            event_mask = ~binary_dilation(
                ~event_mask, iterations=opts.event_mask_buffer_pixels
            )
        logger.debug(f"{int((~event_mask).sum()):,} event-region pixels to fill")

    if opts.unwrap_error_correction:
        disp = correct_region_offset(
            input_disp=disp, mask=mask, wavelength=wavelength_m
        )

    disp = np.ma.filled(disp, np.nan)
    gnss_los = np.ma.filled(gnss_los, np.nan)

    disp_for_fit = (
        fill_masked_region(disp, event_mask) if event_mask is not None else disp
    )

    if downsample_factor > 1:
        logger.info(
            f"Downsampling by factor {downsample_factor} using '{downsample_method}'"
        )

        def _down(arr: np.ndarray) -> np.ndarray:
            return downsample_array(
                arr,
                downsample_factor,
                method=downsample_method,
                weights=downsample_weights,
            )

        disp_for_fit = _down(disp_for_fit)
        gnss_los = _down(gnss_los)
        gnss_los_std = _down(gnss_los_std) if gnss_los_std is not None else None
        win_x = max(1, win_x // downsample_factor)
        win_y = max(1, win_y // downsample_factor)

    # Automatic event detection, only when no event mask was supplied. Runs
    # on the downsampled grid, where single noisy pixels no longer trigger it.
    mad_threshold = opts.residual_outlier_mad_threshold
    region_mad_threshold = opts.residual_region_mad_threshold
    n_auto: int | None = None
    n_region: int | None = None
    if event_mask is None and (
        mad_threshold is not None or region_mad_threshold is not None
    ):
        residual = disp_for_fit - gnss_los
        outlier = np.zeros(residual.shape, dtype=bool)
        if mad_threshold is not None:
            median_res = np.nanmedian(residual)
            # 1.4826 * MAD estimates sigma for Gaussian data (Rousseeuw & Croux
            # 1993), so mad_threshold acts as a robust N-sigma cutoff.
            scaled_mad = np.nanmedian(np.abs(residual - median_res)) * 1.4826
            mad_outlier = np.abs(residual - median_res) > mad_threshold * scaled_mad
            n_auto = int(mad_outlier.sum())
            outlier |= mad_outlier
        if region_mad_threshold is not None:
            region_outlier = detect_coherent_residual_regions(
                residual, region_mad_threshold, opts.residual_region_min_pixels
            )
            n_region = int(region_outlier.sum())
            outlier |= region_outlier
        if opts.event_mask_buffer_pixels > 0:
            outlier = binary_dilation(outlier, iterations=opts.event_mask_buffer_pixels)
        if outlier.any():
            disp_for_fit = fill_masked_region(disp_for_fit, (~outlier).astype(np.uint8))
        logger.debug(f"Automatic event detection filled {int(outlier.sum()):,} pixels")

    sigma = opts.calibration_surface_smoothing_sigma
    if sigma is None:
        # 1/8 of the window suppresses seams without over-smoothing.
        smoothing_sigma: float | None = min(win_x, win_y) / 8
    else:
        smoothing_sigma = sigma or None

    if fit_surface is None:
        from .spatial.processor import SpatialProcessor

        fit_surface = SpatialProcessor().fit_windowed_surface

    logger.debug("Fitting calibration surface...")
    with fit_lock or nullcontext():
        surface = fit_surface(
            insar_data=disp_for_fit,
            gnss_los=gnss_los,
            gnss_los_std=gnss_los_std,
            window_size_x=win_x,
            window_size_y=win_y,
            # 50% overlap makes the Hann-tapered windows sum to ~uniform weight.
            window_overlap_x=win_x // 2,
            window_overlap_y=win_y // 2,
            poly_order=1.5,
            n_jobs=n_jobs,
            smoothing_sigma=smoothing_sigma,
            smoothing_method=opts.calibration_surface_smoothing_method,
            sg_window_length=opts.savitzky_golay.window_length,
            sg_polyorder=opts.savitzky_golay.polyorder,
            mask_residual_outliers=opts.mask_fit_residual_outliers,
        )

    if downsample_factor > 1:
        surface = upsample_array(surface, disp.shape)

    # raw_disp - surface == corrected, reference-zeroed disp - fit.
    # Same precision as the input, so saved surfaces keep their file size.
    surface = np.asarray(surface, dtype=disp.dtype)
    surface += total_correction
    surface += disp0
    return CalibrationSurface(surface, n_auto, n_region)
