# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Fit weights: coherence^p x local robust weights (PRD R-S4, plan T31).

The gamma fit discarded the 15 % most extreme residuals at each end of the
*global* distribution before fitting (~30 % of the pixels per window, real
deformation included); on the F08882 golden pair that mask alone cost about
4.5 mm RMSE at the stations. The v0.5 fit keeps every pixel and weights it:

- **coherence**: ``temporal_coherence ** p`` (p = 8 recommended by the trade
  studies, TS-B1 may revise). The GNSS tie is then made where the GNSS
  stations are: coherent, built-up pixels, not fading rural DS;
- **robust**: pixels that disagree with their *local* low-pass level by more
  than a few scaled MADs are down-weighted (Huber) or gated out, so a single
  unwrapping blunder or a local step cannot pull the kernel.

Ported from the trade study ``surface_methods._robust_weights`` and
``surface_v2.loclin_w``.
"""

from __future__ import annotations

import logging
from typing import Literal

import numpy as np
from scipy.ndimage import gaussian_filter

logger = logging.getLogger(__name__)

__all__ = ["RobustMethod", "coherence_weights", "fit_weights", "robust_weights"]

RobustMethod = Literal["huber", "gate"]
MAD_SCALE = 1.4826


def coherence_weights(
    coherence: np.ndarray | None, power: float, shape: tuple[int, int]
) -> np.ndarray:
    """Return ``clip(coherence, 0, 1) ** power`` (all ones for power 0 or no coherence)."""
    if power < 0:
        msg = f"power must be >= 0, got {power}"
        raise ValueError(msg)
    if power == 0 or coherence is None:
        return np.ones(shape, dtype=np.float64)
    coh = np.asarray(coherence, dtype=np.float64)
    if coh.shape != tuple(shape):
        msg = f"coherence {coh.shape} does not match the field shape {tuple(shape)}"
        raise ValueError(msg)
    coh = np.clip(np.nan_to_num(coh, nan=0.0), 0.0, 1.0)
    return coh**power


def _local_level(field: np.ndarray, weights: np.ndarray, sigma_px: float) -> np.ndarray:
    """Gaussian-weighted local mean of `field` (NaN where the kernel sees no weight)."""
    num = gaussian_filter(
        np.where(weights > 0, field, 0.0) * weights, sigma_px, mode="nearest"
    )
    den = gaussian_filter(weights, sigma_px, mode="nearest")
    return np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)


def robust_weights(
    field: np.ndarray,
    weights: np.ndarray,
    sigma_px: float,
    *,
    method: RobustMethod = "huber",
    threshold: float = 4.0,
    huber_c: float = 1.345,
    iterations: int = 2,
) -> np.ndarray:
    """Down-weight pixels that disagree with their local low-pass level.

    Parameters
    ----------
    field : np.ndarray
        Residual being fitted (finite where ``weights > 0``).
    weights : np.ndarray
        Starting weights (base x coherence); zero excludes a pixel.
    sigma_px : float
        Kernel sigma of the fit (`venti.calibration.loclin.kernel_sigma_px`).
    method : {"huber", "gate"}
        ``"huber"``: continuous weights ``min(1, huber_c * MAD / |d|)`` on the
        standardised local residual d; ``"gate"``: 0/1 at ``threshold`` MADs
        (the trade-study behaviour).
    threshold : float
        Gate width in scaled MADs for ``"gate"`` and the hard cut applied on
        top of Huber (a residual beyond it gets weight 0).
    huber_c : float
        Huber tuning constant (1.345 = 95 % efficiency at the normal).
    iterations : int
        Re-estimate the local level with the updated weights this many times.

    Returns
    -------
    np.ndarray
        Multiplicative robust weights in [0, 1] (1 where nothing is suspicious).

    """
    f = np.asarray(field, dtype=np.float64)
    w0 = np.asarray(weights, dtype=np.float64)
    if f.shape != w0.shape:
        msg = f"field {f.shape} and weights {w0.shape} differ"
        raise ValueError(msg)
    if iterations < 1:
        msg = "iterations must be >= 1"
        raise ValueError(msg)
    rw = np.ones_like(w0)
    for _ in range(iterations):
        w = w0 * rw
        level = _local_level(f, w, sigma_px)
        d = f - level
        ok = np.isfinite(d) & (w > 0)
        if ok.sum() < 10:
            break
        med = float(np.median(d[ok]))
        mad = MAD_SCALE * float(np.median(np.abs(d[ok] - med)))
        if mad <= 0:
            break
        z = np.abs(d - med) / mad
        if method == "gate":
            rw = np.where(ok & (z <= threshold), 1.0, 0.0)
        else:
            with np.errstate(divide="ignore", invalid="ignore"):
                huber = np.where(z > huber_c, huber_c / z, 1.0)
            rw = np.where(ok & (z <= threshold), huber, 0.0)
        rw = np.where(w0 > 0, rw, 0.0)
    n_down = int(np.sum((rw < 1.0) & (w0 > 0)))
    logger.debug(
        "robust_weights(%s): %d of %d pixels down-weighted",
        method,
        n_down,
        int((w0 > 0).sum()),
    )
    return rw


def fit_weights(
    field: np.ndarray,
    base: np.ndarray,
    sigma_px: float,
    *,
    coherence: np.ndarray | None = None,
    coherence_power: float = 0.0,
    robust: bool = False,
    robust_method: RobustMethod = "huber",
    robust_threshold: float = 4.0,
    iterations: int = 2,
) -> np.ndarray:
    """Return the combined weights ``base x coherence^p x robust``.

    With ``coherence_power=0`` and ``robust=False`` this returns `base`
    unchanged, which is what the gamma-equivalent defaults in
    `calibration_options.weights` select.
    """
    base = np.asarray(base, dtype=np.float64)
    w = base * coherence_weights(coherence, coherence_power, base.shape)
    if robust:
        w = w * robust_weights(
            field,
            w,
            sigma_px,
            method=robust_method,
            threshold=robust_threshold,
            iterations=iterations,
        )
    return w
