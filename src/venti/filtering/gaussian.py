# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Gaussian smoothing utilities."""

from __future__ import annotations

import numpy as np
from skimage.filters import gaussian


def apply_gaussian(array: np.ndarray, sigma: float) -> np.ndarray:
    """Apply Gaussian smoothing to a 2-D array.

    Parameters
    ----------
    array : np.ndarray
        2-D input array.
    sigma : float
        Standard deviation for the Gaussian kernel in pixels.

    Returns
    -------
    np.ndarray
        Smoothed array.

    """
    return gaussian(array, sigma)
