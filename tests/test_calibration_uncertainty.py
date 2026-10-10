# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""sigma_CAL (PRD R-E1, R-E2; plan T35)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.calibration.loclin import kernel_sigma_px, loclin_surface
from venti.calibration.remove_restore import sigma_inflation_inside
from venti.calibration.two_pass import calibrate_pair
from venti.calibration.uncertainty import (
    effective_n,
    fit_sigma,
    resolve_k,
    sigma_cal,
)
from venti.frames import load_frame_table
from venti.workflow.config import (
    AlgorithmParameters,
    CalibrationOptions,
    UncertaintyOptions,
)

PIXEL_M = 300.0
CUTOFF_M = 10_000.0
SIGMA_PX = kernel_sigma_px(CUTOFF_M, PIXEL_M)  # 6.2 px
CYCLE_M = 0.02773


def interior(shape, margin_px):
    m = np.zeros(shape, bool)
    m[margin_px:-margin_px, margin_px:-margin_px] = True
    return m


class TestEffectiveN:
    def test_unit_weights_give_four_pi_sigma_squared(self):
        w = np.ones((128, 128))
        s0, n_eff = effective_n(w, SIGMA_PX)
        inner = interior(w.shape, int(4 * SIGMA_PX))
        np.testing.assert_allclose(s0[inner], 1.0, atol=1e-3)  # 4-sigma truncation
        np.testing.assert_allclose(n_eff[inner], 4 * np.pi * SIGMA_PX**2, rtol=0.01)
        # at the edge the kernel sees half the pixels
        assert n_eff[0, 64] < 0.6 * n_eff[64, 64]
        with pytest.raises(ValueError, match="positive"):
            effective_n(w, 0)

    def test_zero_weights_reduce_n_eff_nearby(self):
        w = np.ones((128, 128))
        w[:, 64:] = 0.0
        _, n_eff = effective_n(w, SIGMA_PX)
        assert n_eff[64, 20] > 1.8 * n_eff[64, 63]
        assert n_eff[64, 120] == 0.0


class TestFitSigma:
    def test_white_noise_standard_error(self):
        """sigma_fit = sigma_noise / sqrt(n_eff) in the interior, within 10 %."""
        rng = np.random.default_rng(0)
        noise = 0.005
        field = noise * rng.standard_normal((160, 160))
        w = np.ones_like(field)
        surface, _ = loclin_surface(field, w, PIXEL_M, CUTOFF_M)
        sigma, n_eff = fit_sigma(field, w, SIGMA_PX, surface)
        inner = interior(field.shape, int(4 * SIGMA_PX))
        expected = noise / np.sqrt(4 * np.pi * SIGMA_PX**2)
        assert np.median(sigma[inner]) == pytest.approx(expected, rel=0.10)
        assert np.median(n_eff[inner]) == pytest.approx(
            4 * np.pi * SIGMA_PX**2, rel=0.02
        )
        # the surface error itself is of that size
        assert np.std(surface[inner]) == pytest.approx(expected, rel=0.25)

    def test_nan_where_the_kernel_sees_nothing_and_errors(self):
        field = np.zeros((64, 64))
        w = np.ones_like(field)
        w[:, 32:] = 0.0
        sigma, n_eff = fit_sigma(field, w, SIGMA_PX, field)
        assert np.isnan(sigma[32, 60])
        assert np.isnan(n_eff[32, 60])
        assert np.isfinite(sigma[32, 10])
        with pytest.raises(ValueError, match="differ"):
            fit_sigma(field, w[:10], SIGMA_PX, field)
        with pytest.raises(ValueError, match="non-negative"):
            fit_sigma(field, -w, SIGMA_PX, field)


class TestCombination:
    def test_formula_and_broadcasting(self):
        grid = np.full((4, 4), 0.001)
        fit = np.full((4, 4), 0.0005)
        s = sigma_cal(grid, fit, 3.9, sigma_tropo=0.0002, sigma_ref=0.0001)
        expected = np.sqrt((3.9 * 0.001) ** 2 + 0.0005**2 + 0.0002**2 + 0.0001**2)
        np.testing.assert_allclose(s, expected)
        # scalars only need a shape; None terms are zero
        np.testing.assert_allclose(sigma_cal(0.001, None, 1.0, shape=(2, 3)), 0.001)
        grid[0, 0] = np.nan
        assert np.isnan(sigma_cal(grid, fit, 1.0)[0, 0])
        with pytest.raises(ValueError, match="shape"):
            sigma_cal(grid, fit[:2], 1.0)
        with pytest.raises(ValueError, match="k must be"):
            sigma_cal(grid, fit, -1.0)
        with pytest.raises(ValueError, match="at least one"):
            sigma_cal(0.1, 0.1, 1.0)

    def test_inflation_inside_areas(self):
        mask = np.zeros((40, 40), bool)
        mask[10:30, 10:30] = True
        infl = sigma_inflation_inside(mask, 5.0, 3.0)
        s = sigma_cal(np.full(mask.shape, 0.001), None, 1.0, inflation=infl)
        np.testing.assert_allclose(s[~mask], 0.001)
        assert s[20, 20] > 2.0 * 0.001
        with pytest.raises(ValueError, match=">= 1"):
            sigma_cal(
                np.ones(mask.shape), None, 1.0, inflation=np.full(mask.shape, 0.5)
            )

    def test_resolve_k(self):
        assert resolve_k(UncertaintyOptions()) == 1.0
        assert resolve_k(UncertaintyOptions(k_grid=3.9)) == 3.9
        with pytest.raises(ValueError, match="frame-parameter table"):
            resolve_k(UncertaintyOptions(k_grid="frame_table"))
        # the table resolves it
        base = AlgorithmParameters()
        base.calibration_options.uncertainty.k_grid = "frame_table"
        applied = load_frame_table().apply(base, 8882)
        assert resolve_k(applied.calibration_options.uncertainty) == 3.9


def _opts(k=1.0, inflate=True):
    return CalibrationOptions(
        surface={
            "method": "loclin",
            "cutoff_wavelength_meters": CUTOFF_M,
            "two_pass": False,
        },
        weights={"coherence_power": 0.0, "robust": False},
        uncertainty={"k_grid": k, "inflate_inside_areas": inflate},
    )


class TestCalibrated:
    def test_z_scores_are_calibrated(self):
        """T35.3: std of (surface - truth) / sigma_CAL in [0.8, 1.25].

        Each realisation has white DISP noise (the fit term) and a smooth GNSS
        grid error whose true sigma is k x the formal sigma (the grid term);
        both terms are of the same size so the test checks both.
        """
        rng = np.random.default_rng(7)
        ny = nx = 96
        k = 3.9
        formal = 0.0001  # 0.1 mm formal grid sigma -> 0.39 mm true
        noise = 0.004  # 4 mm DISP noise -> fit term 4 / sqrt(4 pi 6.2^2) = 0.18 mm
        yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
        inner = interior((ny, nx), int(3 * SIGMA_PX))
        ref = (ny // 2, nx // 2)
        z_all = []
        for _ in range(40):
            truth = 0.003 * xx / nx - 0.002 * yy / ny  # the calibration to find
            # grid error: a random plane whose per-pixel std is k*formal
            # (u, v uniform on [-1, 1] have variance 1/3)
            a, b, c = k * formal / np.sqrt(3) * rng.standard_normal(3)
            u = (xx - nx / 2) / (nx / 2)
            v = (yy - ny / 2) / (ny / 2)
            grid_err = a + b * u * np.sqrt(3) + c * v * np.sqrt(3)
            gnss = -grid_err  # the GNSS field is wrong by grid_err
            disp = truth + noise * rng.standard_normal((ny, nx))
            res = calibrate_pair(
                disp,
                gnss,
                np.ones((ny, nx), bool),
                ref,
                _opts(k=k),
                PIXEL_M,
                CYCLE_M,
                gnss_los_std=np.full((ny, nx), formal),
            )
            # the surface is referenced; the grid error passes into it in
            # full (a plane is reproduced exactly), which is what the k term
            # accounts for
            est = res.cal_gnss_surface + res.cal_reference_offset
            err = est - truth
            z_all.append((err / res.sigma_cal)[inner].ravel())
        z = np.concatenate(z_all)
        assert np.isfinite(z).all()
        assert 0.8 < np.std(z) < 1.25, np.std(z)
        assert abs(np.mean(z)) < 0.15

    def test_result_fields_and_inflation(self):
        rng = np.random.default_rng(1)
        ny = nx = 96
        disp = 0.004 * rng.standard_normal((ny, nx))
        gnss = np.zeros((ny, nx))
        exclude = np.zeros((ny, nx), bool)
        exclude[30:66, 30:66] = True
        std = np.full((ny, nx), 0.0003)
        res = calibrate_pair(
            disp,
            gnss,
            np.ones((ny, nx), bool),
            (10, 10),
            _opts(k=2.0),
            PIXEL_M,
            CYCLE_M,
            gnss_los_std=std,
            exclude_mask=exclude,
            sigma_tropo=0.0002,
        )
        assert res.sigma_cal is not None
        assert res.sigma_cal.shape == disp.shape
        assert res.sigma_cal.dtype == disp.dtype
        assert res.sigma_fit is not None
        assert res.n_eff is not None
        outside = ~exclude
        floor = np.sqrt((2.0 * 0.0003) ** 2 + 0.0002**2)
        assert np.nanmin(res.sigma_cal[outside]) >= floor - 1e-9
        # grows inside the interpolated area
        assert np.nanmedian(res.sigma_cal[exclude]) > 1.3 * np.nanmedian(
            res.sigma_cal[outside]
        )
        flat = calibrate_pair(
            disp,
            gnss,
            np.ones((ny, nx), bool),
            (10, 10),
            _opts(k=2.0, inflate=False),
            PIXEL_M,
            CYCLE_M,
            gnss_los_std=std,
            exclude_mask=exclude,
        )
        assert np.nanmedian(flat.sigma_cal[exclude]) < 1.3 * np.nanmedian(
            flat.sigma_cal[outside]
        )

    def test_without_grid_sigma_only_the_fit_term(self, caplog):
        rng = np.random.default_rng(2)
        disp = 0.004 * rng.standard_normal((64, 64))
        with caplog.at_level("WARNING", logger="venti.calibration.two_pass"):
            res = calibrate_pair(
                disp,
                np.zeros_like(disp),
                np.ones(disp.shape, bool),
                (5, 5),
                _opts(),
                PIXEL_M,
                CYCLE_M,
            )
        assert "fit term only" in caplog.text
        np.testing.assert_allclose(res.sigma_cal, res.sigma_fit, rtol=1e-6)

    def test_frame_table_k_must_be_resolved(self):
        o = _opts()
        o.uncertainty.k_grid = "frame_table"
        disp = np.zeros((32, 32))
        with pytest.raises(ValueError, match="frame-parameter table"):
            calibrate_pair(
                disp, disp, np.ones(disp.shape, bool), (3, 3), o, PIXEL_M, CYCLE_M
            )

    def test_gamma_path_has_no_sigma(self):
        disp = np.zeros((40, 40))
        o = CalibrationOptions(
            unwrap_error_correction=False, window_size_meters=20 * PIXEL_M
        )
        res = calibrate_pair(
            disp, disp, np.ones(disp.shape, bool), (3, 3), o, PIXEL_M, CYCLE_M, n_jobs=1
        )
        assert res.sigma_cal is None
        assert res.sigma_fit is None


def test_sigma_maps_cover_a_frame_that_is_not_a_multiple_of_the_factor():
    """Regression: a 97 x 101 frame at factor 6 (the F08882 frame is 7733 x 9464)."""
    rng = np.random.default_rng(3)
    disp = 0.004 * rng.standard_normal((97, 101))
    o = _opts()
    o.downsample_factor = 6
    res = calibrate_pair(
        disp,
        np.zeros_like(disp),
        np.ones(disp.shape, bool),
        (10, 10),
        o,
        PIXEL_M,
        CYCLE_M,
        gnss_los_std=np.full(disp.shape, 0.0003),
    )
    assert res.sigma_cal is not None
    assert res.sigma_cal.shape == disp.shape
    assert res.sigma_fit is not None
    assert res.sigma_fit.shape == disp.shape
    assert res.n_eff is not None
    assert res.n_eff.shape == disp.shape
    assert np.isfinite(res.sigma_cal).all()
    # the trailing rows/columns beyond the last full block take the nearest node
    np.testing.assert_array_equal(res.sigma_fit[96], res.sigma_fit[95])
    np.testing.assert_array_equal(res.sigma_fit[:, 100], res.sigma_fit[:, 99])


def test_float32_frames_give_float32_maps():
    """Full-frame maps stay in the displacement's precision (memory, T37)."""
    rng = np.random.default_rng(4)
    disp = (0.004 * rng.standard_normal((97, 101))).astype(np.float32)
    o = _opts()
    o.downsample_factor = 3
    res = calibrate_pair(
        disp,
        np.zeros_like(disp),
        np.ones(disp.shape, bool),
        (10, 10),
        o,
        PIXEL_M,
        CYCLE_M,
        gnss_los_std=np.full(disp.shape, 0.0003, dtype=np.float32),
    )
    for name in ("sigma_cal", "sigma_fit", "n_eff", "coverage", "cal_gnss_surface"):
        assert getattr(res, name).dtype == np.float32, name
    # float64 callers keep float64
    assert sigma_cal(np.full((2, 2), 0.001), None, 1.0).dtype == np.float64
    assert (
        sigma_cal(np.full((2, 2), 0.001, np.float32), 0.0005, 2.0).dtype == np.float32
    )


def test_without_diagnostic_maps_the_product_layers_are_identical():
    """`diagnostic_maps=False` drops coverage / sigma_fit / n_eff (memory, T47)
    and changes neither the calibration nor sigma_CAL."""
    rng = np.random.default_rng(5)
    disp = (0.004 * rng.standard_normal((97, 101))).astype(np.float32)
    o = _opts()
    o.downsample_factor = 6
    args = (
        disp,
        np.zeros_like(disp),
        np.ones(disp.shape, bool),
        (10, 10),
        o,
        PIXEL_M,
        CYCLE_M,
    )
    kw = {"gnss_los_std": np.full(disp.shape, 0.0003, np.float32)}
    full = calibrate_pair(*args, **kw)
    lean = calibrate_pair(*args, **kw, diagnostic_maps=False)
    assert lean.coverage is None
    assert lean.sigma_fit is None
    assert lean.n_eff is None
    np.testing.assert_array_equal(lean.calibration, full.calibration)
    np.testing.assert_array_equal(lean.sigma_cal, full.sigma_cal)
    # the reference offset is one constant, held as a read-only view
    ref = lean.cal_reference_offset
    assert ref.shape == disp.shape
    assert not ref.flags.writeable
    assert np.all(ref == lean.reference_value)
    lean.assert_closed()
