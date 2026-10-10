# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Weighted least-squares solver (port of MATLAB ``lscov``)."""

from __future__ import annotations

import logging

import numpy as np
import scipy.linalg
import scipy.sparse as spu

logger = logging.getLogger(__name__)


def _weighted_lscov(
    design_matrix: np.ndarray,
    b: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Weighted least-squares via QR factorisation (port of MATLAB ``lscov``).

    Parameters
    ----------
    design_matrix : np.ndarray
        Design matrix of shape ``(n_obs, n_params)``.
    b : np.ndarray
        Observation vector of length ``n_obs``.
    weights : np.ndarray, optional
        Per-observation weights (inverse-variance). Uniform by default.

    Returns
    -------
    x : np.ndarray
        Estimated parameters, shape ``(n_params, 1)``.
    stdx : np.ndarray
        Standard errors of the parameters.
    mse : float
        Mean squared error.
    Qxx : np.ndarray
        Covariance matrix of the parameters.
    res : np.ndarray
        Unweighted residuals.
    wres : np.ndarray
        Weighted residuals.
    obs_hat : np.ndarray
        Fitted observations.
    Qcov : np.ndarray
        Aposteriori cofactor matrix.

    """
    if weights is None:
        weights = np.ones(b.shape)

    b = b[:, np.newaxis]
    n_obs, n_x = design_matrix.shape
    n_r = n_obs - n_x

    if spu.issparse(design_matrix):
        design_matrix = design_matrix.toarray()

    Aw = design_matrix * np.sqrt(weights[:, np.newaxis])
    Bw = b * np.sqrt(weights[:, np.newaxis])

    Q, R, perm = scipy.linalg.qr(Aw, mode="economic", pivoting=True)
    z = Q.T @ Bw

    r_diag = np.diag(R)
    keep = np.abs(r_diag) > np.abs(r_diag[0]) * max(n_obs, n_x) * np.finfo(R.dtype).eps
    rank = int(keep.sum())
    if rank < n_x:
        logger.warning("Design matrix rank-deficient: %d / %d columns kept", rank, n_x)
        R = R[np.ix_(keep, keep)]
        z = z[keep, :]
        perm = perm[keep]

    xx = np.linalg.lstsq(R, z, rcond=None)[0]
    x = np.zeros((n_x, 1))
    x[perm] = xx

    Q_mat = Q if rank == n_x else Q[:, keep]
    wres = Bw - Q_mat @ z
    mse = float(np.sum(wres * wres.conj()) / n_r) if n_r > 0 else 0.0

    Rinv = np.triu(np.linalg.lstsq(R, np.eye(rank), rcond=None)[0])
    Qxx = np.zeros((n_x, n_x))
    Qxx[np.ix_(perm, perm)] = Rinv @ Rinv.T

    stdx = np.sqrt(mse * np.diag(Qxx))
    res = wres / np.sqrt(weights[:, np.newaxis])
    obs_hat = Q_mat @ z / np.sqrt(weights[:, np.newaxis])
    Qcov = (
        design_matrix @ (Qxx / mse) @ design_matrix.T
        if n_r > 0
        else np.zeros(Qxx.shape)
    )

    return x, stdx, mse, Qxx, res, wres, obs_hat, Qcov
