"""SpatialProcessor: high-level interface for InSAR calibration surface fitting."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from ..filtering.moving_window import fit_windowed_plane

logger = logging.getLogger(__name__)


@dataclass
class SpatialProcessor:
    """Wrapper around windowed polynomial plane fitting for InSAR calibration.

    Provides a single :meth:`fit_windowed_surface` method that estimates the
    long-wavelength bias between an InSAR displacement field and a GNSS LOS
    reference surface using a moving-window least-squares polynomial fit.

    Examples
    --------
    ::

        processor = SpatialProcessor()
        surface = processor.fit_windowed_surface(
            insar_data=disp_mm,
            gnss_los=gnss_los_mm,
            window_size_x=1000,
            window_size_y=1000,
            window_overlap_x=10,
            window_overlap_y=10,
            poly_order=1.5,
        )
        calibrated = disp_mm - surface

    """

    def fit_windowed_surface(
        self,
        insar_data: np.ndarray,
        gnss_los: np.ndarray,
        window_size_x: int,
        window_size_y: int,
        window_overlap_x: int = 10,
        window_overlap_y: int = 10,
        window_extend_x: int | None = None,
        window_extend_y: int | None = None,
        gnss_los_std: np.ndarray | None = None,
        poly_order: float = 1.5,
        n_jobs: int = -1,
        smoothing_sigma: float | None = None,
        smoothing_method: str = "gaussian",
        sg_window_length: int = 51,
        sg_polyorder: int = 3,
        mask_residual_outliers: bool = True,
    ) -> np.ndarray:
        """Estimate the long-wavelength InSAR calibration surface.

        Fits a windowed polynomial plane to the difference between
        `insar_data` and `gnss_los`, returning the estimated bias
        surface. Apply the correction as ``calibrated = insar_data - surface``.

        Parameters
        ----------
        insar_data : np.ndarray
            2-D InSAR displacement or velocity array.
        gnss_los : np.ndarray
            2-D GNSS LOS reference field (same shape as `insar_data`).
        window_size_x : int
            Window width in pixels.
        window_size_y : int
            Window height in pixels.
        window_overlap_x : int, optional
            Window overlap in the x direction, by default ``10``.
        window_overlap_y : int, optional
            Window overlap in the y direction, by default ``10``.
        window_extend_x : int, optional
            Extra columns of context added around each window for fitting.
            Defaults to `window_size_x` when not given.
        window_extend_y : int, optional
            Extra rows of context added around each window for fitting.
            Defaults to `window_size_y` when not given.
        gnss_los_std : np.ndarray, optional
            Uncertainty of the GNSS LOS field, same shape as `gnss_los`.
            When given, used as inverse-variance weights in the windowed
            least-squares fit. ``None`` (default) fits with uniform weights.
        poly_order : float, optional
            Polynomial order for plane fitting, by default ``1.5``.
        n_jobs : int, optional
            Number of parallel ``joblib`` workers, by default ``-1`` (all CPUs).
        smoothing_sigma : float, optional
            Standard deviation (pixels) passed to the post-assembly low-pass
            filter.  Ignored when ``smoothing_method="savitzky_golay"``.
            ``None`` disables smoothing.
        smoothing_method : str, optional
            Post-assembly low-pass filter.  One of ``"gaussian"`` (default),
            ``"gaussian_fft"``, ``"hanning_fft"``, or ``"savitzky_golay"``.
        sg_window_length : int, optional
            Window length for Savitzky-Golay in pixels (must be odd),
            by default ``51``.
        sg_polyorder : int, optional
            Polynomial order for Savitzky-Golay, by default ``3``.
        mask_residual_outliers : bool, optional
            Whether `fit_windowed_plane` excludes the most extreme 15% / 85%
            of ``insar_data - gnss_los`` residual pixels before filling and
            fitting, by default ``True``. This is a fixed quantile trim,
            unconditional and independent of any caller-level outlier
            detection — see `fit_windowed_plane`'s docstring. ``False``
            fits `insar_data` as given, masking only genuinely invalid
            (NaN) pixels.

        Returns
        -------
        np.ndarray
            Estimated calibration surface to subtract from `insar_data`,
            shape ``(ny, nx)``.

        """
        win_extend_x = window_extend_x if window_extend_x is not None else window_size_x
        win_extend_y = window_extend_y if window_extend_y is not None else window_size_y

        surface, _ = fit_windowed_plane(
            insar_data=insar_data,
            gnss_los=gnss_los,
            win_xsize=window_size_x,
            win_ysize=window_size_y,
            win_overlap_x=window_overlap_x,
            win_overlap_y=window_overlap_y,
            win_extend_x=win_extend_x,
            win_extend_y=win_extend_y,
            gnss_los_std=gnss_los_std,
            poly_order=poly_order,
            n_jobs=n_jobs,
            smoothing_sigma=smoothing_sigma,
            smoothing_method=smoothing_method,
            sg_window_length=sg_window_length,
            sg_polyorder=sg_polyorder,
            mask_residual_outliers=mask_residual_outliers,
        )
        return surface
