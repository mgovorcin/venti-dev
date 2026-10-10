# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Two-pass calibration and component bookkeeping (PRD R-S1, D6; plan T33)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.calibration.two_pass import COMPONENTS, CalibrationResult, calibrate_pair
from venti.surface import estimate_calibration_surface
from venti.workflow.config import CalibrationOptions

PIXEL_M = 180.0  # m (6x downsampled DISP posting, used at factor 1 here)
CYCLE_M = 0.02773  # Sentinel-1 lambda/2, metres
NY, NX = 160, 200


def scene(seed=0, island=True, bowl=True):
    """Synthetic pair in metres: GNSS plane + long-wavelength ramp + bowl + 1-cycle island + noise."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:NY, 0:NX].astype(float)
    x_m, y_m = xx * PIXEL_M, yy * PIXEL_M
    gnss = 0.004 * x_m / (NX * PIXEL_M) - 0.002 * y_m / (
        NY * PIXEL_M
    )  # 4 mm / 2 mm across
    ramp = 0.010 * (x_m + y_m) / ((NX + NY) * PIXEL_M) + 0.003  # the error to calibrate
    bowl_mask = np.zeros((NY, NX), bool)
    bowl_field = np.zeros((NY, NX))
    if bowl:
        r2 = (x_m - 60 * PIXEL_M) ** 2 + (y_m - 60 * PIXEL_M) ** 2
        bowl_field = -0.15 * np.exp(
            -r2 / (2 * (6 * PIXEL_M) ** 2)
        )  # 15 cm bowl, 1 km sigma
        bowl_mask[40:80, 40:80] = True
    island_mask = np.zeros((NY, NX), bool)
    if island:
        island_mask[110:150, 140:190] = True
    valid = np.ones((NY, NX), bool)
    # water all around the island block, so it is a region of its own
    valid[100:160, 125:200] = False
    valid[island_mask] = True
    valid[:, :8] = False  # open sea
    disp = gnss + ramp + bowl_field + 0.0005 * rng.standard_normal((NY, NX))
    disp = disp + np.where(island_mask, CYCLE_M, 0.0)  # +1 cycle unwrapping error
    disp[~valid] = np.nan
    coh = np.clip(rng.uniform(0.5, 0.95, (NY, NX)), 0, 1)
    return disp, gnss, valid, ramp, bowl_mask, island_mask, coh


def opts_v05(two_pass=True, unwrap=False, coh_power=2.0):
    return CalibrationOptions(
        unwrap_error_correction=unwrap,
        surface={
            "method": "loclin",
            "cutoff_wavelength_meters": 12_000.0,
            "fill_gaps": True,
            "two_pass": two_pass,
        },
        weights={"coherence_power": coh_power, "robust": True},
    )


REF = (20, 100)


class TestComponents:
    def test_calibration_is_the_sum_of_components_and_closed(self):
        disp, gnss, valid, _ramp, bowl_mask, _, coh = scene(island=False)
        tropo = np.full(disp.shape, 0.002)
        set_ = np.full(disp.shape, -0.0007)
        res = calibrate_pair(
            disp + tropo + set_,
            gnss,
            valid,
            REF,
            opts_v05(),
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            tropo=tropo,
            set_correction=set_,
            exclude_mask=bowl_mask,
        )
        total = sum(getattr(res, c) for c in COMPONENTS)
        np.testing.assert_array_equal(res.calibration, total)
        res.assert_closed()
        np.testing.assert_array_equal(res.cal_tropo, tropo)
        np.testing.assert_array_equal(res.cal_set, set_)
        assert res.components_applied == {
            "cal_gnss_surface": True,
            "cal_reference_offset": True,
            "cal_tropo": True,
            "cal_set": True,
            "cal_unwrap_shift": False,
        }
        assert set(res.components()) == set(COMPONENTS)
        assert np.all(res.cal_unwrap_shift == 0)
        assert res.n_excluded_pixels == bowl_mask.sum()
        assert res.passes == 2

    def test_reference_offset_is_the_referenced_pixel_and_not_applied_twice(self):
        disp, gnss, valid, _ramp, _, _, coh = scene(island=False, bowl=False)
        res = calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            opts_v05(two_pass=False),
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
        )
        assert res.reference_value == pytest.approx(float(disp[REF]))
        assert np.all(res.cal_reference_offset == res.reference_value)
        calibrated = disp - res.calibration
        # the calibrated field is tied to GNSS: no residual offset at the reference
        assert abs(float((calibrated - gnss)[REF])) < 0.001

    def test_invalid_reference_raises_and_shapes_are_checked(self):
        disp, gnss, valid, *_ = scene()
        with pytest.raises(ValueError, match="Reference pixel"):
            calibrate_pair(
                disp, gnss, valid, (120, 130), opts_v05(), PIXEL_M, CYCLE_M
            )  # in the channel
        with pytest.raises(ValueError, match="shape"):
            calibrate_pair(disp, gnss[:10], valid, REF, opts_v05(), PIXEL_M, CYCLE_M)


class TestLoclinEndToEnd:
    def test_recovers_the_ramp_outside_the_bowl_with_a_mocked_unwrap_hook(self):
        """Plane + bowl + 1-cycle island + noise: with the true shift the surface error is < 1 mm RMS."""
        disp, gnss, valid, _ramp, bowl_mask, island_mask, coh = scene()
        calls = {}

        def hook(residual, fit_valid, cycle_m):
            # a perfect detector: the island is one cycle high
            calls["cycle_m"] = cycle_m
            calls["shape"] = residual.shape
            return np.where(island_mask, cycle_m, 0.0), {"regions": 1}

        res = calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            opts_v05(two_pass=True, unwrap=True),
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            exclude_mask=bowl_mask,
            unwrap_hook=hook,
        )
        assert calls["cycle_m"] == CYCLE_M
        assert res.components_applied["cal_unwrap_shift"]
        assert res.unwrap_decisions == {"regions": 1}
        np.testing.assert_allclose(res.cal_unwrap_shift[island_mask], CYCLE_M)
        assert np.all(res.cal_unwrap_shift[~island_mask] == 0)
        # DISP - calibration should equal the GNSS field plus the bowl plus noise
        calibrated = disp - res.calibration
        err = (calibrated - gnss)[valid & ~bowl_mask]
        assert np.sqrt(np.mean(err**2)) < 0.001, np.sqrt(np.mean(err**2))
        # the bowl survives in the product (remove-restore)
        assert np.nanmin((calibrated - gnss)[bowl_mask]) < -0.10
        assert res.passes == 2
        assert res.fit_residual_std is not None
        assert res.fit_residual_std < 0.001
        assert res.coverage is not None
        assert res.coverage.shape == disp.shape
        res.assert_closed()

    def test_without_unwrap_the_island_offset_pulls_the_surface(self):
        """Why R-U1 exists: an island wider than the kernel cannot be told from signal."""
        disp, gnss, valid, _ramp, bowl_mask, island_mask, coh = scene()
        o = opts_v05(two_pass=True, unwrap=False)
        res = calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            o,
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            exclude_mask=bowl_mask,
        )
        assert not res.components_applied["cal_unwrap_shift"]
        hooked = calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            opts_v05(two_pass=True, unwrap=True),
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            exclude_mask=bowl_mask,
            unwrap_hook=lambda _r, _v, c: (np.where(island_mask, c, 0.0), None),
        )
        # the surface is pulled by the uncorrected cycle at the island ...
        pulled = np.abs(res.cal_gnss_surface - hooked.cal_gnss_surface)
        assert pulled[island_mask].max() > 0.005
        # ... while far from it both agree and the product is still tied to GNSS
        far = valid & ~bowl_mask
        far[80:, 100:] = False
        assert pulled[far].max() < 0.001
        err_far = ((disp - res.calibration) - gnss)[far]
        assert np.sqrt(np.mean(err_far**2)) < 0.0015

    def test_single_pass_loclin_and_downsampling(self):
        disp, gnss, valid, _ramp, bowl_mask, _, coh = scene(island=False)
        o = opts_v05(two_pass=False)
        o.downsample_factor = 2
        res = calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            o,
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            exclude_mask=bowl_mask,
        )
        assert res.passes == 1
        assert res.cal_gnss_surface.shape == disp.shape
        calibrated = disp - res.calibration
        err = (calibrated - gnss)[valid & ~bowl_mask]
        assert np.sqrt(np.mean(err**2)) < 0.0015


class TestGammaPath:
    def test_windowed_plane_reproduces_estimate_calibration_surface(self):
        """method='windowed_plane' with unwrap off equals the gamma estimator exactly."""
        disp, gnss, valid, *_ = scene(island=False, bowl=False)
        o = CalibrationOptions(
            unwrap_error_correction=False, window_size_meters=60 * PIXEL_M
        )
        assert o.surface.method == "windowed_plane"
        tropo = np.full(disp.shape, 0.0015)
        res = calibrate_pair(
            disp + tropo, gnss, valid, REF, o, PIXEL_M, CYCLE_M, tropo=tropo, n_jobs=1
        )
        gamma = estimate_calibration_surface(
            disp + tropo,
            gnss,
            valid & np.isfinite(disp),
            REF,
            60,
            corrections=[tropo],
            options=o,
            wavelength_m=2 * CYCLE_M,
            n_jobs=1,
        )
        np.testing.assert_allclose(res.calibration, gamma.surface, atol=1e-7)
        assert res.method == "windowed_plane"
        assert res.passes == 1
        assert res.coverage is None
        res.assert_closed()

    def test_windowed_plane_with_unwrap_keeps_the_shift_as_a_component(self):
        disp, gnss, valid, _ramp, _, island_mask, _ = scene(bowl=False)
        o = CalibrationOptions(
            unwrap_error_correction=True, window_size_meters=60 * PIXEL_M
        )
        res = calibrate_pair(disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M, n_jobs=1)
        # the gamma corrector segments on the mask: the island (cut off by the
        # channel) is its own region and is shifted by one cycle
        assert res.components_applied["cal_unwrap_shift"]
        assert np.median(res.cal_unwrap_shift[island_mask]) == pytest.approx(
            CYCLE_M, abs=1e-6
        )
        res.assert_closed()


def test_result_dataclass_round_trip():
    z = np.zeros((2, 2))
    r = CalibrationResult(
        cal_gnss_surface=z + 1,
        cal_reference_offset=z + 2,
        cal_tropo=z,
        cal_set=z,
        cal_unwrap_shift=z + 0.5,
        components_applied=dict.fromkeys(COMPONENTS, True),
        method="loclin",
        reference_point=(0, 0),
        reference_value=2.0,
    )
    np.testing.assert_array_equal(r.calibration, 3.5)
    r.assert_closed()


class TestFidelityFixes:
    def test_gamma_path_uses_gnss_sigma_only_when_asked(self):
        """cal-disp passes gnss_los_std to the gamma fit only with weight_fit_by_gnss_uncertainty."""
        disp, gnss, valid, *_ = scene(island=False, bowl=False)
        std = 0.001 + 0.002 * np.linspace(0, 1, disp.shape[1])[None, :].repeat(
            disp.shape[0], 0
        )
        o = CalibrationOptions(
            unwrap_error_correction=False, window_size_meters=60 * PIXEL_M
        )
        plain = calibrate_pair(disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M, n_jobs=1)
        with_std = calibrate_pair(
            disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M, gnss_los_std=std, n_jobs=1
        )
        np.testing.assert_array_equal(plain.calibration, with_std.calibration)
        o.weight_fit_by_gnss_uncertainty = True
        weighted = calibrate_pair(
            disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M, gnss_los_std=std, n_jobs=1
        )
        assert np.abs(weighted.calibration - plain.calibration).max() > 1e-6

    def test_coarse_pass_one_tie_matches_the_fine_one(self, monkeypatch):
        """The near-planar tie computed on a coarse grid agrees with the full-grid one."""
        from venti.calibration import two_pass as tp

        disp, gnss, valid, _ramp, bowl_mask, _island_mask, coh = scene()
        seen = {}

        def hook(residual, fit_valid, cycle_m):
            seen["residual"] = residual.copy()
            return np.zeros_like(residual), None

        o = opts_v05(two_pass=True, unwrap=True)
        monkeypatch.setattr(tp, "TIE_MAX_PX", 10_000)
        calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            o,
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            exclude_mask=bowl_mask,
            unwrap_hook=hook,
        )
        fine = seen["residual"]
        monkeypatch.setattr(tp, "TIE_MAX_PX", 50)  # 160 x 200 -> factor 4
        calibrate_pair(
            disp,
            gnss,
            valid,
            REF,
            o,
            PIXEL_M,
            CYCLE_M,
            coherence=coh,
            exclude_mask=bowl_mask,
            unwrap_hook=hook,
        )
        coarse = seen["residual"]
        ok = valid & ~bowl_mask
        assert np.nanmax(np.abs(fine - coarse)[ok]) < 0.0005


@pytest.mark.parametrize("method", ["windowed_plane", "loclin"])
def test_inputs_are_not_modified(method):
    """calibrate_pair no longer copies its inputs defensively (memory, T37);
    it must never write into them."""
    disp, gnss, valid, _ramp, bowl_mask, _, coh = scene(island=False)
    disp = disp.astype(np.float32)
    gnss32 = gnss.astype(np.float32)
    tropo = np.full(disp.shape, 0.002, dtype=np.float32)
    set_ = np.full(disp.shape, -0.0007, dtype=np.float32)
    before = [a.copy() for a in (disp, gnss32, tropo, set_)]
    if method == "loclin":
        o = opts_v05(two_pass=True)
    else:
        o = CalibrationOptions(
            unwrap_error_correction=False, window_size_meters=60 * PIXEL_M
        )
    res = calibrate_pair(
        disp,
        gnss32,
        valid,
        REF,
        o,
        PIXEL_M,
        CYCLE_M,
        coherence=coh,
        tropo=tropo,
        set_correction=set_,
        exclude_mask=bowl_mask,
        n_jobs=1,
    )
    for a, b in zip((disp, gnss32, tropo, set_), before, strict=True):
        np.testing.assert_array_equal(a, b)
    res.assert_closed()
    # components that are not applied share one read-only zero layer
    assert not res.cal_unwrap_shift.flags.writeable
