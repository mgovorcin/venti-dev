# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Post-assembly low-pass smoothing filters for calibration surfaces.

Each function takes a pre-assembled calibration surface and a boolean validity
mask and returns a smoothed surface of the same shape.  The boundary-aware
pattern (divide smoothed values by smoothed weights) ensures that valid pixels
near the edge of the data footprint are not pulled toward zero.
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)


def gaussian_spatial(
    surface: np.ndarray,
    valid: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Apply a spatial-domain Gaussian low-pass filter to a calibration surface.

    Parameters
    ----------
    surface : np.ndarray
        2-D calibration surface array.
    valid : np.ndarray
        Boolean mask of pixels that contain valid data.
    sigma : float
        Gaussian standard deviation in pixels.

    Returns
    -------
    np.ndarray
        Filtered surface, same shape and dtype as `surface`.

    """
    from scipy.ndimage import gaussian_filter

    _kw: dict = {"sigma": sigma, "truncate": 2.0}
    valid_f = valid.astype(np.float64)
    smoothed_weight = gaussian_filter(valid_f, **_kw)
    has_weight = smoothed_weight > 0
    smoothed_vals = gaussian_filter(surface * valid_f, **_kw)
    result = surface.copy()
    result[has_weight] = smoothed_vals[has_weight] / smoothed_weight[has_weight]
    logger.debug("gaussian_spatial: sigma=%.1f px", sigma)
    return result


def gaussian_fft(
    surface: np.ndarray,
    valid: np.ndarray,
    sigma: float,
) -> np.ndarray:
    r"""Apply a frequency-domain Gaussian low-pass filter to a calibration surface.

    Applies a Gaussian low-pass kernel in the Fourier domain.  Equivalent to
    `gaussian_spatial` but avoids truncation of the Gaussian kernel, so it is
    more accurate for large sigma values at the cost of being a global
    (non-local) operation.

    Parameters
    ----------
    surface : np.ndarray
        2-D calibration surface array.
    valid : np.ndarray
        Boolean mask of pixels that contain valid data.
    sigma : float
        Gaussian standard deviation in pixels (same units as `gaussian_spatial`).

    Returns
    -------
    np.ndarray
        Filtered surface, same shape and dtype as `surface`.

    Notes
    -----
    The frequency-domain transfer function is:

    .. math::

        H(f) = \\exp\\!\\left(-2\\pi^2 \\sigma^2 f^2\\right)

    where :math:`f` is in cycles per pixel.

    """
    from scipy.fft import fft2, fftfreq, ifft2

    ny, nx = surface.shape
    fy = fftfreq(ny)[:, np.newaxis]
    fx = fftfreq(nx)[np.newaxis, :]
    kernel = np.exp(-2.0 * np.pi**2 * sigma**2 * (fy**2 + fx**2))

    valid_f = valid.astype(np.float64)
    smoothed_vals = np.real(ifft2(fft2(surface * valid_f) * kernel))
    smoothed_weight = np.real(ifft2(fft2(valid_f) * kernel))

    has_weight = smoothed_weight > 0
    result = surface.copy()
    result[has_weight] = smoothed_vals[has_weight] / smoothed_weight[has_weight]
    logger.debug("gaussian_fft: sigma=%.1f px", sigma)
    return result


def hanning_fft(
    surface: np.ndarray,
    valid: np.ndarray,
    sigma: float,
) -> np.ndarray:
    r"""Apply a frequency-domain Hanning low-pass filter to a calibration surface.

    Applies a raised-cosine (Hanning) low-pass window in the Fourier domain.
    The Hanning taper has a sharper roll-off than a Gaussian and produces no
    ringing, making it well-suited for surfaces that contain sharp spatial
    gradients near the edges of valid data.

    Parameters
    ----------
    surface : np.ndarray
        2-D calibration surface array.
    valid : np.ndarray
        Boolean mask of pixels that contain valid data.
    sigma : float
        Cutoff half-width in pixels.  The -3 dB frequency is
        ``f_c = 1 / (2 * sigma)`` cycles per pixel.

    Returns
    -------
    np.ndarray
        Filtered surface, same shape and dtype as `surface`.

    Notes
    -----
    The frequency-domain transfer function is:

    .. math::

        H(f) = \\begin{cases}
            \\frac{1}{2}\\left(1 + \\cos\\!\\left(\\pi f / f_c\\right)\\right)
            & |f| \\leq f_c \\\\
            0 & |f| > f_c
        \\end{cases}

    where :math:`f_c = 1 / (2\\sigma)`.

    """
    from scipy.fft import fft2, fftfreq, ifft2

    ny, nx = surface.shape
    fy = fftfreq(ny)[:, np.newaxis]
    fx = fftfreq(nx)[np.newaxis, :]
    f_mag = np.sqrt(fy**2 + fx**2)
    f_c = 1.0 / (2.0 * sigma)
    kernel = np.where(
        f_mag <= f_c,
        0.5 * (1.0 + np.cos(np.pi * f_mag / f_c)),
        0.0,
    )

    valid_f = valid.astype(np.float64)
    smoothed_vals = np.real(ifft2(fft2(surface * valid_f) * kernel))
    smoothed_weight = np.real(ifft2(fft2(valid_f) * kernel))

    has_weight = smoothed_weight > 0
    result = surface.copy()
    result[has_weight] = smoothed_vals[has_weight] / smoothed_weight[has_weight]
    logger.debug("hanning_fft: sigma=%.1f px (f_c=%.4f cy/px)", sigma, f_c)
    return result


def savitzky_golay(
    surface: np.ndarray,
    valid: np.ndarray,
    window_length: int,
    polyorder: int,
) -> np.ndarray:
    """Apply a 2-D Savitzky-Golay low-pass filter to a calibration surface.

    Applies `scipy.signal.savgol_filter` separably along rows then columns.
    Compared to Gaussian smoothing, Savitzky-Golay better preserves local
    curvature and peaks in the surface because it fits a local polynomial
    rather than computing a weighted mean.

    Parameters
    ----------
    surface : np.ndarray
        2-D calibration surface array.
    valid : np.ndarray
        Boolean mask of pixels that contain valid data.
    window_length : int
        Length of the filter window in pixels (must be odd and greater than
        `polyorder`).
    polyorder : int
        Order of the polynomial used to fit the samples within each window.

    Returns
    -------
    np.ndarray
        Filtered surface, same shape and dtype as `surface`.

    """
    from scipy.signal import savgol_filter

    assert window_length % 2 == 1, "`window_length` must be odd"
    assert polyorder < window_length, "`polyorder` must be less than `window_length`"

    # Replace invalid pixels with zero for the filter pass, then renormalise
    # using the same separable weight approach.
    valid_f = valid.astype(np.float64)
    masked = surface * valid_f

    smoothed_vals = savgol_filter(
        savgol_filter(masked, window_length, polyorder, axis=0),
        window_length,
        polyorder,
        axis=1,
    )
    smoothed_weight = savgol_filter(
        savgol_filter(valid_f, window_length, polyorder, axis=0),
        window_length,
        polyorder,
        axis=1,
    )

    has_weight = smoothed_weight > 0
    result = surface.copy()
    result[has_weight] = smoothed_vals[has_weight] / smoothed_weight[has_weight]
    logger.debug(
        "savitzky_golay: window_length=%d, polyorder=%d", window_length, polyorder
    )
    return result
