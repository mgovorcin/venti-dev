# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""The calibration uncertainty sigma_CAL (PRD R-E1, R-E2; plan T35).

::

    sigma_CAL**2 = (k * sigma_grid)**2 + sigma_fit**2 + sigma_tropo**2 + sigma_ref**2

- ``sigma_grid``: the GNSS LOS sigma of the pair, from the gridded E/N/U
  sigmas (`venti.gnss.sampling.project_field_to_los`). The UNR grid's formal
  sigma is optimistic: on F08882 the grid-minus-MIDAS RMS was 1.75 mm against
  a formal 0.44 mm, hence the inflation ``k`` (3.9 there; per frame from the
  frame-parameter table once TS-G1 fits it; 1 keeps the gamma meaning);
- ``sigma_fit``: how well the kernel determines the surface from the scatter
  of the residual around it: the kernel-weighted residual variance divided
  by the effective number of pixels the kernel sees (`fit_sigma`);
- ``sigma_tropo``, ``sigma_ref``: caller-supplied terms for the tropospheric
  correction and the reference offset (0 until a model exists);
- inside interpolated defo/event areas sigma grows with the distance into
  the area (`venti.calibration.remove_restore.sigma_inflation_inside`);
- unwrap shifts are treated as exact (R-E1).

sigma_CAL covers the *calibration* only (R-E2): DISP's own noise below the
cutoff (~10 mm) is the user's to add. `fit_sigma` assumes the residual
scatter is uncorrelated between pixels; correlated short-wavelength
deformation left in the residual makes it optimistic by the square root of
the correlation area in pixels. All functions keep the units of their inputs.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.ndimage import gaussian_filter

from ..workflow.config import UncertaintyOptions

logger = logging.getLogger(__name__)

__all__ = ["effective_n", "fit_sigma", "resolve_k", "sigma_cal"]


def effective_n(weights: np.ndarray, sigma_px: float) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(s0, n_eff)``: kernel coverage and effective pixel count.

    For a normalised Gaussian kernel ``K`` and pixel weights ``w`` the
    kernel-weighted mean has variance ``sigma**2 * sum(K**2 w**2) /
    sum(K w)**2``, so the effective number of independent pixels is
    ``n_eff = (sum K w)**2 / sum (K**2 w**2)``. ``K**2`` is a Gaussian with
    sigma ``sigma_px / sqrt(2)`` scaled by ``1 / (4 pi sigma_px**2)``, which
    gives ``n_eff = 4 pi sigma_px**2 * s0**2 / G_{sigma/sqrt2}(w**2)``: with
    unit weights everywhere ``n_eff = 4 pi sigma_px**2``.
    """
    w = np.asarray(weights, dtype=np.float64)
    if sigma_px <= 0:
        msg = "sigma_px must be positive"
        raise ValueError(msg)
    s0 = gaussian_filter(w, sigma_px, mode="constant", cval=0.0)
    s2 = gaussian_filter(w * w, sigma_px / np.sqrt(2.0), mode="constant", cval=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        n_eff = np.where(s2 > 0, 4.0 * np.pi * sigma_px**2 * s0 * s0 / s2, 0.0)
    return s0, n_eff


def fit_sigma(
    field: np.ndarray,
    weights: np.ndarray,
    sigma_px: float,
    surface: np.ndarray,
    *,
    min_coverage: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(sigma_fit, n_eff)`` for a local-linear surface.

    Parameters
    ----------
    field : np.ndarray
        The residual that was fitted (gap-filled, finite where ``weights > 0``).
    weights : np.ndarray
        The fit weights (`venti.calibration.weights.fit_weights` output).
    sigma_px : float
        Kernel sigma in pixels (`venti.calibration.loclin.kernel_sigma_px`).
    surface : np.ndarray
        The fitted surface (`venti.calibration.loclin.loclin_surface`).
    min_coverage : float
        Below this kernel coverage the surface is an extrapolation and
        ``sigma_fit`` is NaN (the caller masks or fills it).

    Returns
    -------
    sigma_fit, n_eff : np.ndarray
        Standard error of the surface at each pixel
        (``sqrt(local residual variance / n_eff)``), and the effective
        number of pixels behind it.

    """
    f = np.asarray(field, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    s = np.asarray(surface, dtype=np.float64)
    if not (f.shape == w.shape == s.shape):
        msg = f"field {f.shape}, weights {w.shape} and surface {s.shape} differ"
        raise ValueError(msg)
    if np.any(w < 0):
        msg = "weights must be non-negative"
        raise ValueError(msg)
    r2 = np.where(w > 0, (f - s) ** 2, 0.0)
    s0, n_eff = effective_n(w, sigma_px)
    sr2 = gaussian_filter(w * r2, sigma_px, mode="constant", cval=0.0)
    ok = (s0 > min_coverage) & (n_eff > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        var_local = np.where(ok, sr2 / s0, np.nan)
        sigma = np.sqrt(np.where(ok, var_local / n_eff, np.nan))
    logger.debug(
        "fit_sigma: sigma_px %.1f, n_eff median %.0f, sigma_fit median %.3g",
        sigma_px,
        float(np.nanmedian(np.where(ok, n_eff, np.nan))),
        float(np.nanmedian(sigma)),
    )
    return sigma, np.where(ok, n_eff, np.nan)


def resolve_k(options: UncertaintyOptions) -> float:
    """Return the numeric ``k``; ``'frame_table'`` must have been applied already."""
    k = options.k_grid
    if k == "frame_table":
        msg = (
            "uncertainty.k_grid is 'frame_table': apply the frame-parameter table "
            "(FrameParameterTable.apply) before calibrating"
        )
        raise ValueError(msg)
    if not isinstance(k, (int, float)) or k < 0:
        msg = f"k_grid must be a non-negative number, got {k!r}"
        raise ValueError(msg)
    return float(k)


def sigma_cal(
    sigma_grid: np.ndarray | float | None,
    sigma_fit: np.ndarray | float | None,
    k: float,
    *,
    sigma_tropo: np.ndarray | float | None = None,
    sigma_ref: np.ndarray | float | None = None,
    inflation: np.ndarray | None = None,
    shape: tuple[int, int] | None = None,
) -> np.ndarray:
    """Combine the R-E1 terms into sigma_CAL (units of the inputs).

    Parameters
    ----------
    sigma_grid : array or float, optional
        GNSS LOS sigma of the pair (NaN where the grid has no value).
    sigma_fit : array or float, optional
        Surface standard error (`fit_sigma`).
    k : float
        Inflation of `sigma_grid` (`resolve_k`).
    sigma_tropo, sigma_ref : array or float, optional
        Extra terms, default 0.
    inflation : array, optional
        Multiplicative factor >= 1 inside interpolated areas
        (`sigma_inflation_inside`).
    shape : (ny, nx), optional
        Output shape when every term is a scalar.

    Returns
    -------
    np.ndarray
        ``inflation * sqrt((k sigma_grid)**2 + sigma_fit**2 + sigma_tropo**2
        + sigma_ref**2)``; NaN where any term is NaN.

    """
    if k < 0:
        msg = f"k must be >= 0, got {k}"
        raise ValueError(msg)
    terms = [
        k * np.asarray(0.0 if sigma_grid is None else sigma_grid, dtype=np.float64),
        np.asarray(0.0 if sigma_fit is None else sigma_fit, dtype=np.float64),
        np.asarray(0.0 if sigma_tropo is None else sigma_tropo, dtype=np.float64),
        np.asarray(0.0 if sigma_ref is None else sigma_ref, dtype=np.float64),
    ]
    arrays = [t for t in terms if t.ndim > 0] + (
        [np.asarray(inflation, dtype=np.float64)] if inflation is not None else []
    )
    if arrays:
        out_shape = arrays[0].shape
        for a in arrays:
            if a.shape != out_shape:
                msg = f"sigma terms differ in shape: {[a.shape for a in arrays]}"
                raise ValueError(msg)
    elif shape is not None:
        out_shape = tuple(shape)
    else:
        msg = "give at least one array term or `shape`"
        raise ValueError(msg)
    total = np.zeros(out_shape, dtype=np.float64)
    for t in terms:
        total = total + t * t
    sigma = np.sqrt(total)
    if inflation is not None:
        infl = np.asarray(inflation, dtype=np.float64)
        if np.nanmin(infl) < 1.0:
            msg = "inflation must be >= 1 everywhere"
            raise ValueError(msg)
        sigma = sigma * infl
    return sigma
