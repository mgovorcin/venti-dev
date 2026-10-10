# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Local-linear calibration surface with a physical cutoff (PRD R-S2, plan T30).

The gamma surface is a moving-window plane whose size is not the scale it
resolves (half-response is about 4.5 x the window). Here the surface is a
weighted local-linear regression with a Gaussian kernel whose width follows
from a **half-response wavelength**: the Gaussian smoother's transfer
function is ``exp(-2 pi^2 sigma^2 / lambda^2)``, so

    sigma = lambda_cutoff * sqrt(ln 2 / 2) / pi

attenuates a sinusoid of wavelength ``lambda_cutoff`` by exactly one half.
For uniform coverage the local-linear fit equals the Gaussian-weighted mean
(the slope terms cancel by symmetry), so the cutoff carries over; near
edges and in uneven coverage the linear terms remove the bias a plain
smoother would have.

The regression is solved for every pixel at once from Gaussian-filtered
moments (`local_linear_surface`), which is O(N log N) regardless of the
kernel width. Ported from the trade study
``surface_fitting/continuous_surface/improved/surface_v2.py``.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.ndimage import gaussian_filter

logger = logging.getLogger(__name__)

__all__ = ["kernel_sigma_px", "local_linear_surface", "loclin_surface", "smoothstep"]

# exp(-2 pi^2 sigma^2 / lambda^2) = 1/2  ->  sigma = lambda * sqrt(ln2 / 2) / pi
_HALF_RESPONSE = float(np.sqrt(np.log(2.0) / 2.0) / np.pi)


def kernel_sigma_px(cutoff_wavelength_m: float, pixel_m: float) -> float:
    """Gaussian sigma (pixels) with half response at `cutoff_wavelength_m`."""
    if cutoff_wavelength_m <= 0 or pixel_m <= 0:
        msg = "cutoff_wavelength_m and pixel_m must be positive"
        raise ValueError(msg)
    return cutoff_wavelength_m * _HALF_RESPONSE / pixel_m


def smoothstep(x: np.ndarray) -> np.ndarray:
    """Hermite ramp 0 -> 1 over [0, 1] (3x^2 - 2x^3), clipped outside."""
    t = np.clip(np.asarray(x, dtype=np.float64), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def local_linear_surface(
    field: np.ndarray,
    weights: np.ndarray,
    sigma_px: float,
    ridge: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Weighted local-linear fit of `field` at every pixel with a Gaussian kernel.

    Parameters
    ----------
    field : np.ndarray
        2-D values to fit (finite wherever ``weights > 0``).
    weights : np.ndarray
        Non-negative per-pixel weights (0 excludes a pixel).
    sigma_px : float
        Kernel sigma in pixels (`kernel_sigma_px`).
    ridge : float
        Relative Tikhonov term that keeps the 3x3 normal matrix invertible
        where the kernel sees almost no data.

    Returns
    -------
    surface, coverage : np.ndarray
        The fitted level at each pixel, and the kernel-weighted coverage
        ``S0`` (Gaussian-filtered weights; ~1 inside dense data, -> 0 far
        from any data). Where ``coverage`` is tiny the surface is a
        regularised extrapolation and callers should blend or mask it.

    """
    f = np.asarray(field, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if f.shape != w.shape or f.ndim != 2:
        msg = f"field {f.shape} and weights {w.shape} must be the same 2-D shape"
        raise ValueError(msg)
    if np.any(w < 0):
        msg = "weights must be non-negative"
        raise ValueError(msg)
    if not np.isfinite(f[w > 0]).all():
        msg = "field must be finite wherever weights > 0 (fill gaps first)"
        raise ValueError(msg)

    ny, nx = f.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
    xx /= sigma_px
    yy /= sigma_px
    r = np.where(w > 0, f, 0.0)

    def g(a: np.ndarray) -> np.ndarray:
        return gaussian_filter(a, sigma_px, mode="constant", cval=0.0)

    # Gaussian-weighted moments; the local centring (x - x0) is applied
    # algebraically so one set of filters serves every pixel.
    s0, sx, sy = g(w), g(w * xx), g(w * yy)
    sxx, sxy, syy = g(w * xx * xx), g(w * xx * yy), g(w * yy * yy)
    r0, rx, ry = g(w * r), g(w * xx * r), g(w * yy * r)
    cx, cy = sx - xx * s0, sy - yy * s0
    cxx = sxx - 2 * xx * sx + xx * xx * s0
    cyy = syy - 2 * yy * sy + yy * yy * s0
    cxy = sxy - xx * sy - yy * sx + xx * yy * s0
    normal = np.stack(
        [
            np.stack([s0, cx, cy], -1),
            np.stack([cx, cxx, cxy], -1),
            np.stack([cy, cxy, cyy], -1),
        ],
        -2,
    )
    rhs = np.stack([r0, rx - xx * r0, ry - yy * r0], -1)
    normal = normal + ridge * np.maximum(s0, ridge)[..., None, None] * np.eye(3)
    level = np.linalg.solve(normal, rhs[..., None])[..., 0, 0]
    return level, s0


def loclin_surface(
    field: np.ndarray,
    weights: np.ndarray,
    pixel_m: float,
    cutoff_wavelength_m: float,
    *,
    fallback: np.ndarray | None = None,
    coverage_blend: tuple[float, float] = (0.05, 0.3),
) -> tuple[np.ndarray, np.ndarray]:
    """Calibration surface with a half-response cutoff; continuous everywhere.

    Parameters
    ----------
    field : np.ndarray
        Residual to model (DISP minus GNSS LOS), gap-filled so it is finite.
    weights : np.ndarray
        Per-pixel weights: `venti.calibration.gaps.base_weights` times the
        coherence/robust weights of `venti.calibration.weights`.
    pixel_m : float
        Posting of `field` in meters.
    cutoff_wavelength_m : float
        Half-response wavelength (R-S2; 50 km in the trade studies).
    fallback : np.ndarray, optional
        Field used where the kernel sees almost no weight (open sea,
        outside the swath): by default the Gaussian-smoothed `field` itself.
    coverage_blend : (low, high)
        Coverage range over which the local-linear estimate is blended into
        the fallback with a smoothstep, so the surface stays continuous
        without letting empty regions feed back onto the data.

    Returns
    -------
    surface, coverage : np.ndarray
        Float64 surface (finite everywhere) and the coverage map.

    """
    sigma = kernel_sigma_px(cutoff_wavelength_m, pixel_m)
    level, coverage = local_linear_surface(field, weights, sigma)
    if fallback is None:
        fallback = gaussian_filter(
            np.asarray(field, dtype=np.float64), sigma, mode="nearest"
        )
    low, high = coverage_blend
    beta = smoothstep((coverage - low) / (high - low))
    surface = beta * np.where(beta > 0, level, 0.0) + (1.0 - beta) * fallback
    if not np.isfinite(surface).all():  # pragma: no cover - defensive
        msg = "loclin_surface produced non-finite values"
        raise RuntimeError(msg)
    logger.debug(
        "loclin_surface: cutoff %.0f m -> sigma %.1f px; coverage p5/p50 = %.2f/%.2f",
        cutoff_wavelength_m,
        sigma,
        float(np.percentile(coverage, 5)),
        float(np.percentile(coverage, 50)),
    )
    return surface, coverage
