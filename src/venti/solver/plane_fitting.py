"""Polynomial surface evaluation and uncertainty propagation."""

from __future__ import annotations

import numpy as np

from .design_matrix import _POLY_N_COEFF, _design_matrix_poly
from .lscov import _weighted_lscov


def _calc_plane_values(
    x: np.ndarray,
    y: np.ndarray,
    coef: np.ndarray,
    poly_order: float = 1.0,
) -> np.ndarray:
    """Evaluate a fitted polynomial surface at given coordinates.

    Parameters
    ----------
    x : np.ndarray
        X coordinates (same shape as output).
    y : np.ndarray
        Y coordinates.
    coef : np.ndarray
        Polynomial coefficients from :func:`_weighted_lscov`.
    poly_order : float, optional
        Polynomial order, by default ``1.0``.

    Returns
    -------
    np.ndarray
        Evaluated surface values.

    """
    if poly_order == 0:
        return np.ones(x.shape) * coef[0]
    if poly_order == 1:
        return x * coef[0] + y * coef[1] + coef[2]
    if poly_order == 1.5:
        return x * coef[0] + y * coef[1] + x * y * coef[2] + coef[3]
    if poly_order == 2:
        return (
            x**2 * coef[0]
            + y**2 * coef[1]
            + x * y * coef[2]
            + x * coef[3]
            + y * coef[4]
            + coef[5]
        )
    if poly_order == 3:
        return (
            x**3 * coef[0]
            + y**3 * coef[1]
            + x**2 * coef[2]
            + y**2 * coef[3]
            + x * y * coef[4]
            + x * coef[5]
            + y * coef[6]
        )
    msg = f"Unsupported poly_order={poly_order}"
    raise ValueError(msg)


def _calc_plane_uncertainty(
    x: np.ndarray,
    y: np.ndarray,
    covariance: np.ndarray,
    poly_order: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate coefficient covariance to surface uncertainty.

    Parameters
    ----------
    x : np.ndarray
        X coordinates.
    y : np.ndarray
        Y coordinates.
    covariance : np.ndarray
        Coefficient covariance matrix from :func:`_weighted_lscov`.
    poly_order : float, optional
        Polynomial order, by default ``1.0``.

    Returns
    -------
    plane_var : np.ndarray
        Propagated variance.
    plane_std : np.ndarray
        Propagated standard deviation.

    """
    e = np.ones(x.shape)
    if poly_order == 0:
        plane_var = np.ones(x.shape) * covariance[0, 0]
    elif poly_order == 1:
        x_var, y_var, z_var = np.diag(covariance)
        cxy, cxz, cyz = covariance[0, 1], covariance[0, 2], covariance[1, 2]
        plane_var = (
            x**2 * x_var
            + y**2 * y_var
            + e**2 * z_var
            + 2 * x * (y * cxy + e * cxz)
            + 2 * y * (e * cyz)
        )
    elif poly_order == 1.5:
        x_var, y_var, xy_var, z_var = np.diag(covariance)
        cxy, cxxy, cxz = covariance[0, 1], covariance[0, 2], covariance[0, 3]
        cyxy, cyz, cxyz = covariance[1, 2], covariance[1, 3], covariance[2, 3]
        plane_var = (
            x**2 * x_var
            + y**2 * y_var
            + (x * y) ** 2 * xy_var
            + e**2 * z_var
            + 2 * x * (y * cxy + x * y * cxxy + e * cxz)
            + 2 * y * (x * y * cyxy + e * cyz)
            + 2 * x * y * (e * cxyz)
        )
    else:
        # Higher orders: diagonal-only approximation
        diag = np.diag(covariance)
        A = _design_matrix_poly(x.ravel(), y.ravel(), poly_order=poly_order)
        plane_var = (A**2 @ diag).reshape(x.shape)

    return plane_var, np.sqrt(np.abs(plane_var))


def _fit_plane(
    data: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray,
    order: float = 1.5,
    decimate: int = 1,
    smooth: bool = False,
    smooth_sigma: float = 5.0,
    data_std: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a polynomial surface to 2-D data at given coordinates.

    Parameters
    ----------
    data : np.ndarray
        2-D data array.
    lons : np.ndarray
        2-D array of x/longitude coordinates (same shape as `data`).
    lats : np.ndarray
        2-D array of y/latitude coordinates (same shape as `data`).
    order : float, optional
        Polynomial order, by default ``1.5``.
    decimate : int, optional
        Subsampling factor for the inversion, by default ``1``.
    smooth : bool, optional
        Apply Gaussian smoothing before fitting, by default ``False``.
    smooth_sigma : float, optional
        Sigma for Gaussian smoothing in pixels, by default ``5.0``.
    data_std : np.ndarray, optional
        Per-pixel uncertainty of `data`, same shape, used as inverse-variance
        weights (``w = 1 / data_std**2``) in the least-squares solve.
        Non-positive or non-finite entries fall back to uniform weight for
        that pixel. ``None`` (default) fits with uniform weights.

    Returns
    -------
    plane : np.ndarray
        Fitted surface, same shape as `data`.
    plane_std : np.ndarray
        Propagated uncertainty of the fit.

    """
    if smooth:
        from ..filtering.gaussian import apply_gaussian

        data = apply_gaussian(data, smooth_sigma)

    data_sub = data[::decimate, ::decimate]
    lons_sub = lons[::decimate, ::decimate]
    lats_sub = lats[::decimate, ::decimate]

    valid = ~np.isnan(data_sub.ravel())
    n_coeff = _POLY_N_COEFF[order]
    if valid.sum() < n_coeff:
        return np.zeros_like(data), np.zeros_like(data)

    A = _design_matrix_poly(
        lons_sub.ravel()[valid], lats_sub.ravel()[valid], poly_order=order
    )
    b = data_sub.ravel()[valid]
    if data_std is None:
        w = np.ones_like(b)
    else:
        std_sub = data_std[::decimate, ::decimate].ravel()[valid]
        w = np.ones_like(b)
        finite_positive = np.isfinite(std_sub) & (std_sub > 0)
        w[finite_positive] = 1.0 / std_sub[finite_positive] ** 2

    _, _, _, Qxx, res, *_ = _weighted_lscov(A, b, w)

    res_flat = res.ravel()
    inliers = (res_flat >= res_flat.mean() - 2 * res_flat.std()) & (
        res_flat <= res_flat.mean() + 2 * res_flat.std()
    )
    if inliers.sum() >= n_coeff:
        x1, _, _, Qxx, *_ = _weighted_lscov(A[inliers], b[inliers], w[inliers])
    else:
        x1, _, _, Qxx, *_ = _weighted_lscov(A, b, w)

    plane = _calc_plane_values(lons, lats, x1.ravel(), order)
    _, plane_std = _calc_plane_uncertainty(lons, lats, Qxx, order)
    return plane, plane_std
