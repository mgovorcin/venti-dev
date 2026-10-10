# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Per-pixel LOS -> (E, U) with N from GNSS, and single-geometry projection.

A LOS displacement is ``los = e E + n N + u U`` with ``(e, n, u)`` the unit
look vector of the DISP-S1-STATIC ``line_of_sight_enu`` layers. Two or more
geometries (ascending + descending) separate east from up once the north
component is fixed from the GNSS grid (D15): per pixel

    los_k - n_k N = e_k E + u_k U,      k = 1..K,

solved by weighted least squares with ``w_k = 1 / (sigma_k^2 + n_k^2
sigma_N^2)``; ``sigma_k`` is the look's displacement sigma (sigma_CAL
combined with the DISP noise, R-E3). With one geometry east must come from
GNSS too:

    U = (los - e E_gnss - n N_gnss) / u.

**Validity of the single-geometry projection (plan T53.4).** The GNSS grid
resolves horizontal motion only at wavelengths above its node spacing
(about 25 km, and 50 km after calibration smoothing). Where the horizontal
motion has shorter wavelengths (fault creep, landslides, a subsiding
basin's horizontal flanks) the projected ``U`` absorbs ``e (E - E_gnss) / u``
of the unmodelled east motion, i.e. up to about 0.8 of it for Sentinel-1
(``e ~ 0.6``, ``u ~ 0.75``). The product therefore carries a per-pixel
mode flag (`MODE_PROJECTION`) and documents that projected pixels are
valid for long-wavelength horizontal motion only.

Conditioning: the weighted design's condition number is returned per pixel;
`decompose` falls back to projection where it exceeds `cond_max` (two looks
from the same side are nearly collinear in ``(e, u)``). For a Sentinel-1
ascending + descending pair ``cond ~ 1.2``.

All displacement inputs share one unit (metres in the products); `north`
and `east` are displacements over the pair's span, i.e. GNSS rate x dt.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

__all__ = [
    "MODE_NONE",
    "MODE_PROJECTION",
    "MODE_WLS",
    "DecompositionResult",
    "decompose",
    "decompose_wls",
    "project_vertical",
]

# mode flag values: no valid look / E and U from the WLS solve (>= 2 looks) /
# U from the single-geometry projection with GNSS east (1 look or ill-conditioned)
MODE_NONE = 0
MODE_WLS = 1
MODE_PROJECTION = 2

ENU = tuple[np.ndarray, np.ndarray, np.ndarray]


@dataclass
class DecompositionResult:
    """Per-pixel E/U and their sigmas, the mode and the conditioning."""

    east: np.ndarray
    up: np.ndarray
    sigma_east: np.ndarray
    sigma_up: np.ndarray
    mode: np.ndarray
    cond: np.ndarray
    n_looks: np.ndarray

    def counts(self) -> dict[str, int]:
        return {
            "none": int(np.sum(self.mode == MODE_NONE)),
            "wls": int(np.sum(self.mode == MODE_WLS)),
            "projection": int(np.sum(self.mode == MODE_PROJECTION)),
        }


def _as_field(
    value: np.ndarray | float, shape: tuple[int, ...], name: str
) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim == 0:
        return np.full(shape, float(arr))
    if arr.shape != shape:
        msg = f"{name} shape {arr.shape} does not match the LOS grid {shape}"
        raise ValueError(msg)
    return arr


def _check_looks(
    los: Sequence[np.ndarray],
    sigma: Sequence[np.ndarray | float],
    enu: Sequence[ENU],
    min_looks: int,
) -> tuple[list[np.ndarray], list[np.ndarray], list[ENU], tuple[int, ...]]:
    if len(los) < min_looks:
        msg = f"need at least {min_looks} look(s), got {len(los)}"
        raise ValueError(msg)
    if not (len(los) == len(sigma) == len(enu)):
        msg = "los, sigma and enu must have one entry per look"
        raise ValueError(msg)
    shape = np.shape(los[0])
    los_f = [_as_field(d, shape, f"los[{k}]") for k, d in enumerate(los)]
    sig_f = [_as_field(s, shape, f"sigma[{k}]") for k, s in enumerate(sigma)]
    enu_f: list[ENU] = []
    for k, triple in enumerate(enu):
        if len(triple) != 3:
            msg = f"enu[{k}] must be (east, north, up)"
            raise ValueError(msg)
        e, n, u = (_as_field(c, shape, f"enu[{k}]") for c in triple)
        enu_f.append((e, n, u))
    for k, s in enumerate(sig_f):
        if np.any(s[np.isfinite(s)] < 0):
            msg = f"sigma[{k}] has negative values"
            raise ValueError(msg)
    return los_f, sig_f, enu_f, shape


def decompose_wls(
    los: Sequence[np.ndarray],
    sigma: Sequence[np.ndarray | float],
    enu: Sequence[ENU],
    north: np.ndarray | float,
    *,
    sigma_north: np.ndarray | float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve ``[E, U]`` per pixel from K >= 2 looks with north fixed.

    Parameters
    ----------
    los : sequence of arrays
        LOS displacement per look (NaN where invalid), all on one grid.
    sigma : sequence of arrays or floats
        Displacement sigma per look (sigma_CAL combined with DISP noise).
    enu : sequence of (east, north, up)
        Unit look vectors per look (arrays or floats).
    north : array or float
        North displacement from GNSS over the pair's span.
    sigma_north : array or float, optional
        Its sigma; enters each look's weight as ``n_k^2 sigma_N^2``.

    Returns
    -------
    east, up, sigma_east, sigma_up, cond : np.ndarray
        NaN where fewer than two looks are valid or the system is singular.
        `cond` is the condition number of the weighted design matrix.

    """
    los_f, sig_f, enu_f, shape = _check_looks(los, sigma, enu, min_looks=2)
    n_field = _as_field(north, shape, "north")
    sn_field = _as_field(sigma_north, shape, "sigma_north")

    m11 = np.zeros(shape)
    m12 = np.zeros(shape)
    m22 = np.zeros(shape)
    b1 = np.zeros(shape)
    b2 = np.zeros(shape)
    n_valid = np.zeros(shape, dtype=np.int16)
    for d, s, (e, n, u) in zip(los_f, sig_f, enu_f, strict=True):
        var = s * s + (n * sn_field) ** 2
        ok = np.isfinite(d) & np.isfinite(var) & (var > 0) & np.isfinite(e + n + u)
        w = np.where(ok, 1.0 / np.where(ok, var, 1.0), 0.0)
        rhs = np.where(ok, d - n * n_field, 0.0)
        m11 += w * e * e
        m12 += w * e * u
        m22 += w * u * u
        b1 += w * e * rhs
        b2 += w * u * rhs
        n_valid += ok.astype(np.int16)

    det = m11 * m22 - m12 * m12
    trace = m11 + m22
    with np.errstate(invalid="ignore", divide="ignore"):
        half = 0.5 * trace
        disc = np.sqrt(np.maximum(half * half - det, 0.0))
        lam_max, lam_min = half + disc, half - disc
        cond = np.sqrt(np.where(lam_min > 0, lam_max / lam_min, np.inf))
        # singular when the small eigenvalue is rounding noise of the large one
        solvable = (n_valid >= 2) & np.isfinite(det) & (lam_min > 1e-10 * lam_max)
        inv_det = np.where(solvable, 1.0 / np.where(solvable, det, 1.0), np.nan)
        east = (m22 * b1 - m12 * b2) * inv_det
        up = (m11 * b2 - m12 * b1) * inv_det
        sigma_east = np.sqrt(m22 * inv_det)
        sigma_up = np.sqrt(m11 * inv_det)
    cond = np.where(solvable, cond, np.nan)
    logger.debug(
        "decompose_wls: %d of %d pixels solved, cond median %.2f",
        int(solvable.sum()),
        solvable.size,
        float(np.nanmedian(cond)) if solvable.any() else float("nan"),
    )
    return east, up, sigma_east, sigma_up, cond


def project_vertical(
    los: np.ndarray,
    sigma: np.ndarray | float,
    enu: ENU,
    east: np.ndarray | float,
    north: np.ndarray | float,
    *,
    sigma_east: np.ndarray | float = 0.0,
    sigma_north: np.ndarray | float = 0.0,
    min_up: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Project one geometry to vertical with east and north from GNSS.

    ``U = (los - e E - n N) / u`` and
    ``sigma_U^2 = (sigma^2 + e^2 sigma_E^2 + n^2 sigma_N^2) / u^2``.
    Valid for long-wavelength horizontal motion only (module docstring).
    Pixels with ``|u| < min_up`` are NaN.
    """
    los_f, sig_f, enu_f, shape = _check_looks([los], [sigma], [enu], min_looks=1)
    d, s, (e, n, u) = los_f[0], sig_f[0], enu_f[0]
    e_field = _as_field(east, shape, "east")
    n_field = _as_field(north, shape, "north")
    se = _as_field(sigma_east, shape, "sigma_east")
    sn = _as_field(sigma_north, shape, "sigma_north")
    ok = np.isfinite(d) & (np.abs(u) >= min_up) & np.isfinite(e + n + u)
    with np.errstate(invalid="ignore", divide="ignore"):
        up = np.where(ok, (d - e * e_field - n * n_field) / u, np.nan)
        var = (s * s + (e * se) ** 2 + (n * sn) ** 2) / (u * u)
        sigma_up = np.where(ok, np.sqrt(var), np.nan)
    return up, sigma_up


def decompose(
    los: Sequence[np.ndarray],
    sigma: Sequence[np.ndarray | float],
    enu: Sequence[ENU],
    east_gnss: np.ndarray | float,
    north_gnss: np.ndarray | float,
    *,
    sigma_east_gnss: np.ndarray | float = 0.0,
    sigma_north_gnss: np.ndarray | float = 0.0,
    cond_max: float = 10.0,
    min_up: float = 0.1,
) -> DecompositionResult:
    """Per pixel: WLS where two or more looks are well conditioned, else projection.

    Parameters
    ----------
    los, sigma, enu
        As in `decompose_wls`; one look is allowed.
    east_gnss, north_gnss : array or float
        GNSS east and north displacement over the span (the north is used by
        both modes, the east by the projection only).
    sigma_east_gnss, sigma_north_gnss : array or float, optional
        Their sigmas.
    cond_max : float
        Above this condition number the WLS solution is not trusted and the
        pixel is projected instead (flagged `MODE_PROJECTION`).
    min_up : float
        Minimum ``|u|`` for the projection.

    Returns
    -------
    DecompositionResult
        ``east`` is NaN on projected pixels (not solved); ``mode`` is
        `MODE_NONE`, `MODE_WLS` or `MODE_PROJECTION`; ``cond`` is NaN where
        no WLS solve was possible; ``n_looks`` counts valid looks per pixel.

    """
    los_f, sig_f, enu_f, shape = _check_looks(los, sigma, enu, min_looks=1)
    if cond_max <= 1:
        msg = "cond_max must be > 1"
        raise ValueError(msg)
    n_looks = np.zeros(shape, dtype=np.int16)
    for d in los_f:
        n_looks += np.isfinite(d).astype(np.int16)

    east = np.full(shape, np.nan)
    up = np.full(shape, np.nan)
    sigma_east = np.full(shape, np.nan)
    sigma_up = np.full(shape, np.nan)
    cond = np.full(shape, np.nan)
    mode = np.full(shape, MODE_NONE, dtype=np.uint8)

    if len(los_f) >= 2:
        e_w, u_w, se_w, su_w, c_w = decompose_wls(
            los_f, sig_f, enu_f, north_gnss, sigma_north=sigma_north_gnss
        )
        use = np.isfinite(u_w) & np.isfinite(c_w) & (c_w <= cond_max)
        east[use], up[use] = e_w[use], u_w[use]
        sigma_east[use], sigma_up[use] = se_w[use], su_w[use]
        cond = c_w
        mode[use] = MODE_WLS

    # projection: inverse-variance mean of every valid look's projected up
    num = np.zeros(shape)
    den = np.zeros(shape)
    for d, s, triple in zip(los_f, sig_f, enu_f, strict=True):
        u_k, su_k = project_vertical(
            d,
            s,
            triple,
            east_gnss,
            north_gnss,
            sigma_east=sigma_east_gnss,
            sigma_north=sigma_north_gnss,
            min_up=min_up,
        )
        ok = np.isfinite(u_k) & np.isfinite(su_k) & (su_k > 0)
        w = np.where(ok, 1.0 / np.where(ok, su_k * su_k, 1.0), 0.0)
        num += w * np.where(ok, u_k, 0.0)
        den += w
    proj = (mode == MODE_NONE) & (den > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        up[proj] = (num / den)[proj]
        sigma_up[proj] = np.sqrt(1.0 / den)[proj]
    mode[proj] = MODE_PROJECTION

    result = DecompositionResult(
        east=east,
        up=up,
        sigma_east=sigma_east,
        sigma_up=sigma_up,
        mode=mode,
        cond=cond,
        n_looks=n_looks,
    )
    logger.info("decompose: %s", result.counts())
    return result
