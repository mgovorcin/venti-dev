# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Polynomial design matrix construction for plane fitting."""

from __future__ import annotations

import numpy as np

# Polynomial order -> number of coefficients
_POLY_N_COEFF = {0: 1, 1: 3, 1.5: 4, 2: 6, 3: 8}


def _design_matrix_poly(
    x: np.ndarray,
    y: np.ndarray,
    c: float = 1,
    poly_order: float = 1.0,
) -> np.ndarray:
    """Build the polynomial design matrix for plane fitting.

    Parameters
    ----------
    x : np.ndarray
        1-D array of x (easting/longitude) coordinates.
    y : np.ndarray
        1-D array of y (northing/latitude) coordinates.
    c : float, optional
        Constant offset scale, by default ``1``.
    poly_order : float, optional
        Polynomial order: ``0``, ``1``, ``1.5``, ``2``, or ``3``.

    Returns
    -------
    np.ndarray
        Design matrix of shape ``(n_obs, n_coeff)``.

    """
    ones = np.ones(len(x)) * c
    if poly_order == 0:
        A = np.array([ones])
    elif poly_order == 1:
        A = np.array([x, y, ones])
    elif poly_order == 1.5:
        A = np.array([x, y, x * y, ones])
    elif poly_order == 2:
        A = np.array([x**2, y**2, x * y, x, y, ones])
    elif poly_order == 3:
        A = np.array([x**3, y**3, x**2, y**2, x * y, x, y, ones])
    else:
        msg = f"Unsupported poly_order={poly_order}. Choose from {list(_POLY_N_COEFF)}"
        raise ValueError(msg)
    return A.T
