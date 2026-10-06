# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Raster resampling utilities (downsample and upsample)."""

from __future__ import annotations

import logging
import warnings

import numpy as np
from scipy.ndimage import zoom

logger = logging.getLogger(__name__)


def downsample_array(
    array: np.ndarray,
    factor: int,
    method: str = "mean",
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Downsample array by given factor using specified aggregation method.

    Parameters
    ----------
    array : np.ndarray
        2D array to downsample
    factor : int
        Downsampling factor (e.g., 2 = half resolution)
    method : str, optional
        Aggregation method: 'mean' or 'median', by default 'mean'
    weights : np.ndarray, optional
        Weight array for weighted downsampling (same shape as array).
        If provided, computes weighted mean. Ignored if method='median'.

    Returns
    -------
    np.ndarray
        Downsampled array

    Examples
    --------
    ::

        # Simple mean downsampling
        downsampled = downsample_array(data, factor=4)

        # Median downsampling
        downsampled = downsample_array(data, factor=4, method='median')

        # Weighted mean downsampling
        downsampled = downsample_array(data, factor=4, method='mean', weights=coherence)

    Notes
    -----
    NaN values are handled using nanmean or nanmedian. Blocks with all NaN
    values will result in NaN in the output.

    """
    if factor == 1:
        return array

    if method not in ["mean", "median"]:
        msg = f"Invalid method '{method}'. Must be 'mean' or 'median'"
        raise ValueError(msg)

    new_shape = (array.shape[0] // factor, array.shape[1] // factor)

    trimmed_rows = new_shape[0] * factor
    trimmed_cols = new_shape[1] * factor
    array_trimmed = array[:trimmed_rows, :trimmed_cols]

    if weights is not None and method == "mean":
        weights_trimmed = weights[:trimmed_rows, :trimmed_cols]
        weights_trimmed = np.where(np.isnan(weights_trimmed), 0, weights_trimmed)
        weights_trimmed = np.where(weights_trimmed < 0, 0, weights_trimmed)

    blocks = array_trimmed.reshape(
        new_shape[0], factor, new_shape[1], factor
    ).transpose(0, 2, 1, 3)

    if method == "mean":
        if weights is not None:
            weight_blocks = weights_trimmed.reshape(
                new_shape[0], factor, new_shape[1], factor
            ).transpose(0, 2, 1, 3)

            with np.errstate(invalid="ignore", divide="ignore"):
                weight_blocks_masked = np.where(np.isnan(blocks), 0, weight_blocks)
                data_masked = np.where(np.isnan(blocks), 0, blocks)

                weighted_sum = np.sum(data_masked * weight_blocks_masked, axis=(2, 3))
                weight_sum = np.sum(weight_blocks_masked, axis=(2, 3))

                downsampled = weighted_sum / weight_sum
                downsampled = np.where(weight_sum == 0, np.nan, downsampled)
        else:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", r"Mean of empty slice")
                downsampled = np.nanmean(blocks, axis=(2, 3))
    else:  # median
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
            downsampled = np.nanmedian(blocks, axis=(2, 3))

    logger.debug(
        f"Downsampled array from {array.shape} to {downsampled.shape} "
        f"(factor={factor}, method={method}, weighted={weights is not None})"
    )

    return downsampled


def upsample_array(array: np.ndarray, target_shape: tuple[int, int]) -> np.ndarray:
    """Upsample array to target shape using bilinear interpolation.

    Parameters
    ----------
    array : np.ndarray
        2D array to upsample
    target_shape : tuple
        Target shape (rows, cols)

    Returns
    -------
    np.ndarray
        Upsampled array

    Examples
    --------
    ::

        upsampled = upsample_array(downsampled_data, original_shape)

    Notes
    -----
    NaN values are preserved during upsampling.

    """
    if array.shape == target_shape:
        return array

    zoom_factors = (target_shape[0] / array.shape[0], target_shape[1] / array.shape[1])

    nan_mask = np.isnan(array)
    array_filled = np.where(nan_mask, 0, array)

    upsampled = zoom(array_filled, zoom_factors, order=1)

    mask_upsampled = zoom(nan_mask.astype(float), zoom_factors, order=0) > 0.5
    upsampled[mask_upsampled] = np.nan

    logger.debug(f"Upsampled array from {array.shape} to {upsampled.shape}")

    return upsampled
