# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Local-linear surface with a physical cutoff (PRD R-S2, plan T30)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.calibration.gaps import base_weights, fill_gaps
from venti.calibration.loclin import (
    kernel_sigma_px,
    local_linear_surface,
    loclin_surface,
    smoothstep,
)

PIXEL_M = 180.0  # 6x downsampled DISP posting
CUTOFF_M = 50_000.0


def test_kernel_sigma_from_cutoff():
    sigma = kernel_sigma_px(CUTOFF_M, PIXEL_M)
    # sigma = lambda * sqrt(ln2/2)/pi = 50 km * 0.1874 = 9.37 km = 52.1 px
    assert sigma == pytest.approx(52.05, abs=0.1)
    # the Gaussian transfer function is exactly 1/2 at the cutoff
    assert np.exp(
        -2 * np.pi**2 * (sigma * PIXEL_M) ** 2 / CUTOFF_M**2
    ) == pytest.approx(0.5)
    with pytest.raises(ValueError, match="positive"):
        kernel_sigma_px(0, PIXEL_M)


def _response(wavelength_m: float, n: int = 1024) -> float:
    """Amplitude response of loclin to a sinusoid, measured in the interior."""
    pixel_m = PIXEL_M
    x = np.arange(n) * pixel_m
    field = np.sin(2 * np.pi * x / wavelength_m)[None, :].repeat(64, axis=0)
    w = np.ones_like(field)
    level, _cov = local_linear_surface(field, w, kernel_sigma_px(CUTOFF_M, pixel_m))
    inner = np.s_[20:44, n // 4 : 3 * n // 4]
    return float(
        np.sqrt(np.mean(level[inner] ** 2)) / np.sqrt(np.mean(field[inner] ** 2))
    )


@pytest.mark.parametrize(
    ("wavelength_m", "expected"),
    [
        (CUTOFF_M, 0.5),  # half response at the cutoff, by construction
        (
            200_000.0,
            np.exp(
                -2
                * np.pi**2
                * (CUTOFF_M * np.sqrt(np.log(2) / 2) / np.pi) ** 2
                / 200_000.0**2
            ),
        ),
        (
            20_000.0,
            np.exp(
                -2
                * np.pi**2
                * (CUTOFF_M * np.sqrt(np.log(2) / 2) / np.pi) ** 2
                / 20_000.0**2
            ),
        ),
    ],
)
def test_transfer_function_matches_the_gaussian(wavelength_m, expected):
    assert _response(wavelength_m) == pytest.approx(expected, abs=0.05)


def test_plane_is_reproduced_exactly_and_edges_are_unbiased():
    yy, xx = np.mgrid[0:120, 0:160].astype(float)
    plane = 2.0 + 0.03 * xx - 0.02 * yy
    level, cov = local_linear_surface(
        plane, np.ones_like(plane), kernel_sigma_px(CUTOFF_M, PIXEL_M)
    )
    # a local-*linear* fit is exact on a plane (up to the 1e-6 relative ridge),
    # including at the edges where a plain Gaussian smoother would be biased
    np.testing.assert_allclose(level, plane, atol=1e-3)
    # coverage is the kernel mass that falls on data: highest at the centre,
    # below 1 on a frame only ~2 sigma wide, lower at the corners
    assert cov.max() <= 1.0 + 1e-9
    assert cov[60, 80] > cov[0, 0] > 0
    assert cov[60, 80] == pytest.approx(cov.max(), rel=0.02)


def test_zero_weight_pixels_do_not_steer_the_surface():
    _yy, xx = np.mgrid[0:100, 0:100].astype(float)
    plane = 0.01 * xx
    field = plane.copy()
    w = np.ones_like(field)
    field[40:60, 40:60] = 1e3  # garbage with zero weight
    w[40:60, 40:60] = 0.0
    level, _ = local_linear_surface(field, w, 15.0)
    np.testing.assert_allclose(level, plane, atol=1e-3)


def test_invalid_inputs():
    f = np.ones((10, 10))
    with pytest.raises(ValueError, match="same 2-D shape"):
        local_linear_surface(f, np.ones((5, 5)), 3.0)
    with pytest.raises(ValueError, match="non-negative"):
        local_linear_surface(f, -np.ones((10, 10)), 3.0)
    g = f.copy()
    g[0, 0] = np.nan
    with pytest.raises(ValueError, match="fill gaps first"):
        local_linear_surface(g, np.ones((10, 10)), 3.0)


def test_smoothstep():
    np.testing.assert_allclose(
        smoothstep(np.array([-1, 0, 0.5, 1, 2])), [0, 0, 0.5, 1, 1]
    )


def test_loclin_surface_is_continuous_with_sea_and_gaps():
    """A frame with open sea along one side and an enclosed gap: finite, never 0, blended."""
    rng = np.random.default_rng(0)
    ny, nx = 120, 180
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    truth = 5.0 + 0.02 * xx - 0.01 * yy  # mm, long-wavelength residual
    data = truth + 0.5 * rng.standard_normal((ny, nx))
    valid = np.ones((ny, nx), bool)
    valid[:, :50] = False  # open sea
    valid[50:70, 100:130] = False  # enclosed gap
    data[~valid] = np.nan
    filled, filled_mask = fill_gaps(data, valid)
    w = base_weights(valid, filled_mask)
    surface, cov = loclin_surface(filled, w, PIXEL_M, CUTOFF_M)
    assert np.isfinite(surface).all()
    assert not np.any(surface == 0.0)
    # on land the surface follows the truth to well within the noise
    land = valid.copy()
    land[:, :70] = False  # keep away from the sea blend
    assert np.abs(surface - truth)[land].max() < 0.6
    # the enclosed gap is interpolated through (not 0, close to the plane)
    assert np.abs(surface - truth)[50:70, 100:130].max() < 0.6
    # over the sea the surface stays continuous: no jump at the coast
    coast_jump = np.abs(np.diff(surface[:, 45:55], axis=1)).max()
    assert coast_jump < 0.5
    assert cov[:, :20].mean() < cov[:, 100:].mean()
