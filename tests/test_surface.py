# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tests for `venti.surface.estimate_calibration_surface` (the array-only core)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.surface import estimate_calibration_surface
from venti.workflow.config import CalibrationOptions

NY, NX = 60, 80
MASK = np.ones((NY, NX), bool)
REF = (30, 40)


class _SpyFit:
    """Record the arguments of the windowed fit; return a zero surface."""

    def __call__(self, **kwargs):
        self.kwargs = kwargs
        return np.zeros_like(kwargs["insar_data"], dtype=np.float64)


def _epoch():
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[:NY, :NX]
    gnss = 0.002 * xx / NX
    disp = gnss + 0.004 * yy / NY + 1e-4 * rng.standard_normal((NY, NX))
    return disp.astype("float32"), gnss.astype("float32")


def _run(disp, gnss, spy=None, **kwargs):
    kwargs.setdefault("options", CalibrationOptions(unwrap_error_correction=False))
    return estimate_calibration_surface(
        disp, gnss, MASK, REF, 20, n_jobs=1, fit_surface=spy, **kwargs
    )


def test_removes_ramp_and_keeps_input_precision():
    disp, gnss = _epoch()
    result = _run(disp, gnss)
    assert result.surface.dtype == np.float32
    assert np.nanstd(disp - result.surface - gnss) < 5e-4


def test_corrections_are_removed_before_the_fit_and_added_back():
    disp, gnss = _epoch()
    tropo = np.full_like(disp, 0.01)
    solid_earth_tide = np.full_like(disp, 0.002)
    spy = _SpyFit()

    result = _run(disp, gnss, spy, corrections=[tropo, solid_earth_tide])

    corrected = disp - tropo - solid_earth_tide
    np.testing.assert_allclose(
        spy.kwargs["insar_data"], corrected - corrected[REF], atol=1e-6
    )
    np.testing.assert_allclose(
        result.surface, tropo + solid_earth_tide + corrected[REF], atol=1e-6
    )


def test_downsampling_fits_on_coarse_grid_and_returns_full_grid():
    disp, gnss = _epoch()
    spy = _SpyFit()
    weights = np.random.default_rng(1).uniform(0.2, 1, disp.shape)

    result = _run(disp, gnss, spy, downsample_factor=4, downsample_weights=weights)

    assert spy.kwargs["insar_data"].shape == (NY // 4, NX // 4)
    assert spy.kwargs["window_size_x"] == 20 // 4
    assert result.surface.shape == disp.shape


def test_event_mask_region_and_buffer_are_filled_before_the_fit():
    disp, gnss = _epoch()
    disp[20:30, 30:40] += 0.05  # local deformation
    event_mask = np.ones(disp.shape, bool)
    # The mask misses a 1 px rim; a 2 px buffer (4-connected, so diamond
    # shaped) covers it, corners included.
    event_mask[21:29, 31:39] = False
    spy = _SpyFit()
    options = CalibrationOptions(
        unwrap_error_correction=False, event_mask_buffer_pixels=2
    )

    result = _run(disp, gnss, spy, event_mask=event_mask, options=options)

    fitted = spy.kwargs["insar_data"] + disp[REF]
    # The whole deformation, rim included, is replaced by nearby values.
    assert np.nanmax(np.abs(fitted[20:30, 30:40] - disp[20:30, 30:40])) > 0.04
    assert np.nanmax(fitted[20:30, 30:40] - gnss[20:30, 30:40]) < 0.02
    # Automatic detection does not run when an event mask is given.
    assert result.n_auto_masked_pixels is None


def test_automatic_detection_flags_unmasked_deformation():
    disp, gnss = _epoch()
    disp[20:30, 30:40] += 0.05
    options = CalibrationOptions(
        unwrap_error_correction=False,
        residual_outlier_mad_threshold=5.0,
        residual_region_mad_threshold=3.0,
        residual_region_min_pixels=20,
    )

    result = _run(disp, gnss, _SpyFit(), options=options)

    assert result.n_auto_masked_pixels >= 100
    assert result.n_region_masked_pixels >= 100
    assert _run(disp, gnss, _SpyFit()).n_auto_masked_pixels is None


@pytest.mark.parametrize(("sigma", "expected"), [(None, 20 / 8), (0, None), (1.5, 1.5)])
def test_smoothing_sigma_option(sigma, expected):
    disp, gnss = _epoch()
    spy = _SpyFit()
    options = CalibrationOptions(
        unwrap_error_correction=False, calibration_surface_smoothing_sigma=sigma
    )
    _run(disp, gnss, spy, options=options)
    assert spy.kwargs["smoothing_sigma"] == expected


def test_window_size_accepts_pair_and_numpy_int():
    disp, gnss = _epoch()
    for window in [(20, 10), np.int64(20)]:
        spy = _SpyFit()
        estimate_calibration_surface(
            disp,
            gnss,
            MASK,
            REF,
            window,
            n_jobs=1,
            fit_surface=spy,
            options=CalibrationOptions(unwrap_error_correction=False),
        )
        assert spy.kwargs["window_size_x"] == 20


def test_mismatched_shapes_raise():
    disp, gnss = _epoch()
    with pytest.raises(ValueError, match=r"corrections\[0\] shape"):
        _run(disp, gnss, corrections=[np.zeros((3, 3))])
    with pytest.raises(ValueError, match="gnss_los shape"):
        _run(disp, gnss[:10])


def test_surface_follows_gnss_offset_up_to_the_coast():
    """Regression: masked pixels carried weight in the window blend.

    Their 0s pulled the surface toward 0 near every mask edge, so a GNSS
    reference that differs by a smooth field (e.g. IGS20 vs plate-fixed)
    changed the calibrated land near coasts by less than that field.
    """
    disp, gnss = _epoch()
    land = np.ones(disp.shape, bool)
    land[:, 55:] = False  # sea on the east side
    disp[~land] = np.nan
    offset = 0.009  # like ~9 mm/yr of plate motion in LOS over one year
    options = CalibrationOptions(unwrap_error_correction=False)

    def surface(gnss_los):
        return estimate_calibration_surface(
            disp, gnss_los, land, REF, 20, n_jobs=1, options=options
        ).surface

    change = surface(gnss) - surface(gnss + offset)

    np.testing.assert_allclose(change[land], offset, atol=1e-5)


def test_unwrap_correction_uses_half_wavelength(monkeypatch):
    """Integration: estimate_calibration_surface quantises unwrap offsets in lambda/2."""
    import venti.surface as surface_mod
    from venti import unwrap as unwrap_pkg
    from venti.workflow.config import CalibrationOptions

    seen = {}

    def fake_correct(input_disp, mask, cycle_length, **kwargs):
        seen["cycle_length"] = cycle_length
        return np.asarray(input_disp)

    # surface.py does `from .unwrap import correct_region_offset` inside the
    # call, so the package namespace is where the lookup happens.
    monkeypatch.setattr(unwrap_pkg, "correct_region_offset", fake_correct)
    n = 48
    disp = np.zeros((n, n), dtype=np.float32)
    gnss = np.zeros((n, n), dtype=np.float32)
    mask = np.ones((n, n), dtype=bool)
    surface_mod.estimate_calibration_surface(
        disp,
        gnss,
        mask,
        ref_point=(n // 2, n // 2),
        window_size=n,
        options=CalibrationOptions(unwrap_error_correction=True),
        wavelength_m=0.05546,
        n_jobs=1,
    )
    assert seen["cycle_length"] == pytest.approx(0.05546 / 2)
