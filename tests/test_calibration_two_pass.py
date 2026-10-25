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

    def test_separate_regions_keep_the_island_cycle_off_the_mainland(self):
        """T62: with surface.regions='separate' the island's +1 cycle becomes an
        island offset; the mainland surface matches the island-free scene and
        the island is tied to GNSS."""
        disp, gnss, valid, _ramp, bowl_mask, island_mask, coh = scene()
        clean, *_ = scene(island=False)
        clean[~valid] = np.nan

        def run(d, regions):
            o = opts_v05(two_pass=True, unwrap=False)
            o = o.model_copy(
                update={"surface": o.surface.model_copy(update={"regions": regions})}
            )
            return calibrate_pair(
                d, gnss, valid, REF, o, PIXEL_M, CYCLE_M,
                coherence=coh, exclude_mask=bowl_mask,
            )  # fmt: skip

        sep, ref = run(disp, "separate"), run(clean, "separate")
        none = run(disp, "none")
        mainland_coast = valid & ~bowl_mask & ~island_mask
        mainland_coast[:80, :] = False  # the coast facing the island
        mainland_coast[:, :100] = False
        # the island's cycle does not bend the mainland surface ...
        moved_sep = np.abs(sep.cal_gnss_surface - ref.cal_gnss_surface)
        moved_none = np.abs(none.cal_gnss_surface - ref.cal_gnss_surface)
        assert moved_sep[mainland_coast].max() < 0.001
        assert moved_none[mainland_coast].max() > 3 * moved_sep[mainland_coast].max()
        # ... and the island itself is tied to GNSS through its own offset
        err_island = ((disp - sep.calibration) - gnss)[island_mask]
        assert np.sqrt(np.mean(err_island**2)) < 0.0015
        assert sep.components_applied == none.components_applied
        np.testing.assert_allclose(
            sep.calibration, sum(getattr(sep, c) for c in COMPONENTS), atol=1e-6
        )

    def test_fit_mask_keeps_a_decorrelating_field_out_of_the_fit(self):
        """T62 (Houston/Matagorda): an unstable field carrying a 3 cm bias
        (phase bias over crops) pulls the surface unless surface.fit_mask
        leaves it out; the surface is still evaluated on the field and
        sigma_CAL is not inflated there."""
        disp, gnss, valid, _ramp, bowl_mask, _, coh = scene(island=False)
        field = np.zeros_like(valid)
        field[20:60, 120:170] = True
        disp = disp + np.where(field, 0.03, 0.0)
        clean, *_ = scene(island=False)

        def run(d, mask):
            o = opts_v05(two_pass=True)
            name = None if mask is None else "stable.tif"
            o = o.model_copy(
                update={"surface": o.surface.model_copy(update={"fit_mask": name})}
            )
            return calibrate_pair(
                d, gnss, valid, REF, o, PIXEL_M, CYCLE_M,
                coherence=coh, exclude_mask=bowl_mask, fit_mask=mask,
            )  # fmt: skip

        skip, keep, ref = run(disp, ~field), run(disp, None), run(clean, ~field)
        near = valid & ~bowl_mask & ~field
        near[:, :100] = False
        near[90:, :] = False
        moved_skip = np.abs(skip.cal_gnss_surface - ref.cal_gnss_surface)[near]
        moved_keep = np.abs(keep.cal_gnss_surface - ref.cal_gnss_surface)[near]
        assert moved_skip.max() < 0.002
        assert moved_keep.max() > 3 * moved_skip.max()
        # evaluated on the field (not NaN) and not counted as a deforming area
        assert np.isfinite(skip.cal_gnss_surface[field]).all()
        assert skip.n_excluded_pixels == keep.n_excluded_pixels
        np.testing.assert_array_equal(
            np.isnan(skip.sigma_cal), np.isnan(keep.sigma_cal)
        )

    def test_named_fit_mask_needs_the_array(self):
        disp, gnss, valid, *_ = scene(island=False)
        o = opts_v05()
        o.surface.fit_mask = "stable.tif"
        with pytest.raises(ValueError, match="fit_mask"):
            calibrate_pair(disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M)
        with pytest.raises(ValueError, match="fit_mask shape"):
            calibrate_pair(
                disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M, fit_mask=valid[1:]
            )

    def test_caller_fit_mask_without_a_name_acts_like_a_named_one(self):
        """A mask the caller reads from its own input (cal-disp runconfig file)
        needs no packaged name and fits exactly like the same mask by name."""
        disp, gnss, valid, _ramp, bowl_mask, _, coh = scene(island=False)
        field = np.zeros_like(valid)
        field[20:60, 120:170] = True
        disp = disp + np.where(field, 0.03, 0.0)
        named = opts_v05(two_pass=True)
        named = named.model_copy(
            update={"surface": named.surface.model_copy(update={"fit_mask": "s.tif"})}
        )
        kw = {"coherence": coh, "exclude_mask": bowl_mask, "fit_mask": ~field}
        by_name = calibrate_pair(disp, gnss, valid, REF, named, PIXEL_M, CYCLE_M, **kw)
        by_caller = calibrate_pair(
            disp, gnss, valid, REF, opts_v05(two_pass=True), PIXEL_M, CYCLE_M, **kw
        )
        np.testing.assert_array_equal(by_caller.calibration, by_name.calibration)
        np.testing.assert_array_equal(by_caller.sigma_cal, by_name.sigma_cal)

    def test_island_offset_uses_every_valid_island_pixel_under_a_fit_mask(self):
        """T46 Houston: the fit mask shapes the mainland fit only; an island
        whose pixels are all outside the mask still gets its own offset (from
        every valid pixel), so its cycle stays off the product."""
        disp, gnss, valid, _ramp, bowl_mask, island_mask, coh = scene()
        o = opts_v05(two_pass=True)
        o = o.model_copy(
            update={
                "surface": o.surface.model_copy(
                    update={"regions": "separate", "fit_mask": "stable.tif"}
                )
            }
        )
        res = calibrate_pair(
            disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M,
            coherence=coh, exclude_mask=bowl_mask, fit_mask=~island_mask,
        )  # fmt: skip
        err_island = ((disp - res.calibration) - gnss)[island_mask]
        assert np.sqrt(np.mean(err_island**2)) < 0.0015

    def test_extrapolated_sigma_follows_the_calibration(self):
        """uncertainty.extrapolate_sigma: sigma_CAL has a value wherever the
        calibration has one; in a decorrelated block (no fit weight) it is the
        nearest fitted value, grown with the distance; off: NaN there."""
        disp, gnss, valid, _ramp, bowl_mask, _, coh = scene(island=False)
        coh = coh.copy()
        coh[100:150, 110:170] = 0.02  # decorrelated: coh**8 leaves no weight
        masked = np.zeros_like(valid)
        masked[20:35, 40:60] = True  # DISP masked on land: filled calibration
        valid = valid & ~masked
        disp = np.where(valid, disp, np.nan)
        # the grid sigma has no value on the open-sea strip, where the
        # (gap-filled) calibration does
        grid_std = np.full(disp.shape, 0.001)
        grid_std[:, :5] = np.nan

        def run(extrapolate):
            o = opts_v05(two_pass=True, coh_power=8.0)
            o = o.model_copy(
                update={
                    "uncertainty": o.uncertainty.model_copy(
                        update={"extrapolate_sigma": extrapolate}
                    )
                }
            )
            return calibrate_pair(
                disp, gnss, valid, REF, o, PIXEL_M, CYCLE_M,
                gnss_los_std=grid_std, coherence=coh,
                exclude_mask=bowl_mask,
            )  # fmt: skip

        old, new = run(False), run(True)
        cal = np.isfinite(new.calibration)
        block = np.zeros_like(valid)
        block[115:135, 125:155] = True
        # before: no sigma in the decorrelated block or on the masked land
        assert np.isnan(old.sigma_cal[block]).all()
        assert np.isnan(old.sigma_cal[masked]).all()
        # now: sigma exactly where the calibration is
        np.testing.assert_array_equal(np.isfinite(new.sigma_cal), cal)
        # unchanged where the fit had weight, larger inside the block
        both = np.isfinite(old.sigma_cal)
        np.testing.assert_allclose(new.sigma_cal[both], old.sigma_cal[both])
        assert np.median(new.sigma_cal[block]) > np.median(old.sigma_cal[both])
        np.testing.assert_allclose(new.calibration, old.calibration, equal_nan=True)

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


def test_frame_table_turns_regions_on_for_the_island_frames():
    """T62 edge pass: surface.regions = separate only where the benchmark showed
    island errors (Houston, Los Angeles, Puerto Rico); NYC and others keep none."""
    from venti.frames import load_frame_table
    from venti.workflow.config import AlgorithmParameters

    table = load_frame_table()
    for frame, expected in ((8882, "separate"), (16940, "separate"), (35991, "separate"),
                            (8622, "none"), (8886, "none"), (99, "none")):  # fmt: skip
        applied = table.apply(AlgorithmParameters(), frame)
        assert applied.calibration_options.surface.regions == expected, frame


def test_frame_table_gives_houston_its_packaged_fit_mask():
    """T62 edge pass: only F08882 fits on stable pixels, and its mask ships
    with the package on the frame's DISP grid."""
    import rasterio

    from venti.frames import fit_mask_path, load_frame_table
    from venti.workflow.config import AlgorithmParameters

    table = load_frame_table()
    for frame, expected in ((8882, "F08882_stable_coh07.tif"), (38238, None),
                            (16940, None), (99, None)):  # fmt: skip
        applied = table.apply(AlgorithmParameters(), frame)
        assert applied.calibration_options.surface.fit_mask == expected, frame
    with rasterio.open(fit_mask_path("F08882_stable_coh07.tif")) as src:
        assert (src.count, src.dtypes[0], src.shape) == (1, "uint8", (7733, 9464))
        assert src.crs.to_epsg() == 32615
        kept = src.read(1, out_shape=(773, 946)).mean()
    assert 0.08 < kept < 0.11
    with pytest.raises(FileNotFoundError):
        fit_mask_path("../frame_parameters.json")
    with pytest.raises(FileNotFoundError):
        fit_mask_path("F00000_missing.tif")


class TestRegionPlane:
    """T62 edge pass (Catalina, F16940): a large island far offshore receives the
    mainland surface extrapolated across the water, which can tilt across it; a
    constant offset then fits one end of the island only."""

    @staticmethod
    def _scene(island_shape=(40, 120), seed=3):
        from venti.calibration.two_pass import _region_offsets  # noqa: F401

        rng = np.random.default_rng(seed)
        labels = np.zeros((200, 300), dtype=np.int32)
        labels[:, :150] = 1  # mainland (largest)
        r0, c0 = 80, 170
        h, w = island_shape
        labels[r0 : r0 + h, c0 : c0 + w] = 2
        rows, cols = np.indices(labels.shape)
        # extrapolated surface error on the island: offset + tilt along the island
        truth = np.where(
            labels == 2, 0.010 + 2e-4 * (cols - c0) - 1e-4 * (rows - r0), 0.0
        )
        surface = np.zeros(labels.shape)
        residual = truth + 0.0005 * rng.standard_normal(labels.shape)
        usable = labels > 0
        return labels, surface, residual, usable, truth

    def test_plane_recovers_the_island_tilt_and_offset_does_not(self):
        from venti.calibration.two_pass import _region_offsets

        labels, surface, residual, usable, truth = self._scene()
        island = labels == 2
        off = _region_offsets(surface, residual, labels, usable, 50)
        pla = _region_offsets(
            surface, residual, labels, usable, 50, model="plane", plane_min_pixels=1000
        )
        err_off = np.abs(off - truth)[island]
        err_pla = np.abs(pla - truth)[island]
        assert err_off.max() > 0.01  # the tilt (2.4 cm end to end) is left over
        assert err_pla.max() < 0.0005
        # the mainland is never touched
        np.testing.assert_array_equal(pla[labels == 1], surface[labels == 1])
        np.testing.assert_array_equal(off[labels == 1], surface[labels == 1])

    def test_plane_is_robust_to_outliers(self):
        from venti.calibration.two_pass import _region_offsets

        labels, surface, residual, usable, truth = self._scene()
        island = np.argwhere(labels == 2)
        bad = island[::7]  # ~14% of the island one cycle off
        residual[bad[:, 0], bad[:, 1]] += CYCLE_M
        pla = _region_offsets(
            surface, residual, labels, usable, 50, model="plane", plane_min_pixels=1000
        )
        assert np.abs(pla - truth)[labels == 2].max() < 0.001

    def test_small_region_keeps_the_offset(self):
        from venti.calibration.two_pass import _region_offsets

        labels, surface, residual, usable, _ = self._scene(island_shape=(20, 30))
        off = _region_offsets(surface, residual, labels, usable, 50)
        pla = _region_offsets(
            surface, residual, labels, usable, 50, model="plane", plane_min_pixels=1000
        )
        np.testing.assert_array_equal(pla, off)  # 600 px < 1000: offset only

    def test_default_is_the_offset_model(self):
        surf = CalibrationOptions().surface
        assert surf.region_model == "offset"
        assert surf.region_plane_min_pixels == 1000
