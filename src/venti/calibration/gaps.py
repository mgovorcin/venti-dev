# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Gap filling before the surface fit (PRD R-S3, plan T29).

The gamma surface was 0 on cells the fit never touched (13 % of F08882,
±50 mm), because masked and water pixels carried no value into the windowed
fit. The v0.5 surface fits a *filled* field: every invalid pixel gets a value
borrowed from its neighbourhood, and the fit down-weights those pixels
(`base_weights`, 0.02 by default) so they anchor the kernel without steering
it. The filled values are never written to a product; they only make the
surface continuous.
"""

from __future__ import annotations

import logging
from typing import Literal

import numpy as np
from scipy import ndimage

logger = logging.getLogger(__name__)

__all__ = ["FillMethod", "base_weights", "fill_gaps"]

FillMethod = Literal["nearest", "nearest_smooth", "biharmonic", "idw"]


def _valid_mask(data: np.ndarray, mask: np.ndarray | None) -> np.ndarray:
    valid = np.isfinite(data)
    if mask is not None:
        valid &= np.asarray(mask, dtype=bool)
    return valid


def fill_gaps(
    data: np.ndarray,
    mask: np.ndarray | None = None,
    method: FillMethod = "nearest_smooth",
    smooth_sigma: float = 3.0,
    idw_max_distance: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Fill invalid pixels so the field is finite everywhere.

    Parameters
    ----------
    data : np.ndarray
        2-D field; NaN marks invalid pixels.
    mask : np.ndarray, optional
        Additional validity mask (True = valid), e.g. the recommended mask
        AND the water mask AND NOT the defo/event areas.
    method : {"nearest", "nearest_smooth", "biharmonic", "idw"}
        ``"nearest"`` copies the nearest valid value (fast, blocky);
        ``"nearest_smooth"`` (default) does that and then smooths the filled
        pixels with a Gaussian of `smooth_sigma` pixels, valid pixels kept
        exactly; ``"biharmonic"`` uses scikit-image's inpainting (smooth,
        slower, good for enclosed gaps); ``"idw"`` is GDAL's inverse-distance
        fill via `venti.spatial.gap_filling.fill_gaps`.
    smooth_sigma : float
        Gaussian sigma in pixels for ``"nearest_smooth"``.
    idw_max_distance : int, optional
        Search distance in pixels for ``"idw"`` (GDAL default when None).

    Returns
    -------
    filled, filled_mask : np.ndarray
        The finite field (float64, same shape) and a boolean array that is
        True where a value was filled in.

    Raises
    ------
    ValueError
        If no pixel is valid.

    """
    arr = np.asarray(data, dtype=np.float64)
    if arr.ndim != 2:
        msg = f"fill_gaps expects a 2-D field, got shape {arr.shape}"
        raise ValueError(msg)
    valid = _valid_mask(arr, mask)
    if not valid.any():
        msg = "fill_gaps: no valid pixel to fill from"
        raise ValueError(msg)
    filled_mask = ~valid
    if not filled_mask.any():
        return arr.copy(), filled_mask

    if method in ("nearest", "nearest_smooth"):
        # indices of the nearest valid pixel for every pixel
        idx = ndimage.distance_transform_edt(
            ~valid, return_distances=False, return_indices=True
        )
        filled = arr[tuple(idx)]
        if method == "nearest_smooth" and smooth_sigma > 0:
            # smooth only where we filled; valid pixels keep their values
            smoothed = ndimage.gaussian_filter(filled, smooth_sigma, mode="nearest")
            filled = np.where(filled_mask, smoothed, arr)
    elif method == "biharmonic":
        from skimage.restoration import inpaint_biharmonic

        work = np.where(valid, arr, 0.0)
        filled = inpaint_biharmonic(work, filled_mask)
        filled = np.where(valid, arr, filled)
    elif method == "idw":
        from ..spatial.gap_filling import fill_gaps as _idw_fill

        work = np.where(valid, arr, np.nan)
        filled = np.asarray(
            _idw_fill(work, filled_mask, max_search_distance=idw_max_distance),
            dtype=np.float64,
        )
        filled = np.where(valid, arr, filled)
    else:  # pragma: no cover - Literal prevents it
        msg = f"unknown fill method {method!r}"
        raise ValueError(msg)

    if not np.isfinite(filled).all():
        # a fill method that leaves holes (idw beyond its search distance)
        # gets the nearest-neighbour fallback so the contract holds
        idx = ndimage.distance_transform_edt(
            ~np.isfinite(filled), return_distances=False, return_indices=True
        )
        filled = filled[tuple(idx)]
    logger.debug(
        "fill_gaps(%s): filled %d of %d pixels (%.1f %%)",
        method,
        int(filled_mask.sum()),
        filled_mask.size,
        100.0 * filled_mask.mean(),
    )
    return filled, filled_mask


def base_weights(
    valid: np.ndarray,
    filled_mask: np.ndarray,
    w_filled: float = 0.02,
) -> np.ndarray:
    """Weights of a filled field: 1 where valid, `w_filled` where filled, 0 elsewhere.

    "Elsewhere" is a pixel that is neither valid nor filled, which `fill_gaps`
    never produces; it is kept so callers can zero out regions on purpose.
    """
    valid = np.asarray(valid, dtype=bool)
    filled_mask = np.asarray(filled_mask, dtype=bool)
    if valid.shape != filled_mask.shape:
        msg = f"valid {valid.shape} and filled_mask {filled_mask.shape} differ"
        raise ValueError(msg)
    if not 0.0 <= w_filled <= 1.0:
        msg = f"w_filled must be in [0, 1], got {w_filled}"
        raise ValueError(msg)
    w = np.zeros(valid.shape, dtype=np.float64)
    w[valid] = 1.0
    w[filled_mask & ~valid] = w_filled
    return w
