# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""LOS decomposition and projection (PRD section 3.2, R-E3; plan T53)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.decomposition import (
    MODE_NONE,
    MODE_PROJECTION,
    MODE_WLS,
    decompose,
    decompose_wls,
    project_vertical,
)

NY, NX = 80, 120


def s1_geometry(ascending: bool):
    """Sentinel-1-like unit look vectors with incidence 30-46 deg across range."""
    inc = np.deg2rad(np.linspace(30, 46, NX))[None, :].repeat(NY, 0)
    heading = np.deg2rad(-12.0 if ascending else 192.0)  # flight direction
    # look is to the right of the heading
    look_az = heading + np.pi / 2
    e = np.sin(inc) * np.sin(look_az) * np.ones((NY, NX))
    n = np.sin(inc) * np.cos(look_az) * np.ones((NY, NX))
    u = np.cos(inc) * np.ones((NY, NX))
    sign = -1.0  # displacement toward the satellite is positive
    return sign * e, sign * n, u


def truth():
    yy, xx = np.mgrid[0:NY, 0:NX].astype(float)
    east = 0.004 * xx / NX - 0.001
    north = 0.002 * np.ones((NY, NX))
    up = (
        -0.03 * np.exp(-(((xx - 60) / 15) ** 2 + ((yy - 40) / 12) ** 2))
        + 0.001 * yy / NY
    )
    return east, north, up


def los_of(enu, east, north, up):
    e, n, u = enu
    return e * east + n * north + u * up


class TestWls:
    def test_recovers_east_and_up_exactly_without_noise(self):
        east, north, up = truth()
        asc, desc = s1_geometry(True), s1_geometry(False)
        los = [los_of(asc, east, north, up), los_of(desc, east, north, up)]
        e_hat, u_hat, se, su, cond = decompose_wls(
            los, [0.002, 0.003], [asc, desc], north
        )
        np.testing.assert_allclose(e_hat, east, atol=1e-12)
        np.testing.assert_allclose(u_hat, up, atol=1e-12)
        assert np.all(np.isfinite(cond))
        assert 1.0 < np.median(cond) < 2.0  # asc + desc is well conditioned
        assert np.all(se > 0)
        assert np.all(su > 0)
        # geometry: east is weaker than up for near-polar orbits
        assert np.median(se) > np.median(su)

    def test_sigmas_are_calibrated_with_noise(self):
        rng = np.random.default_rng(0)
        east, north, up = truth()
        asc, desc = s1_geometry(True), s1_geometry(False)
        sig_a, sig_d, sig_n = 0.003, 0.002, 0.0005
        los = [
            los_of(asc, east, north, up) + sig_a * rng.standard_normal((NY, NX)),
            los_of(desc, east, north, up) + sig_d * rng.standard_normal((NY, NX)),
        ]
        north_gnss = north + sig_n * rng.standard_normal((NY, NX))
        e_hat, u_hat, se, su, _ = decompose_wls(
            los, [sig_a, sig_d], [asc, desc], north_gnss, sigma_north=sig_n
        )
        z_e = (e_hat - east) / se
        z_u = (u_hat - up) / su
        assert 0.93 < np.std(z_e) < 1.07
        assert 0.93 < np.std(z_u) < 1.07
        assert abs(np.mean(z_e)) < 0.05
        assert abs(np.mean(z_u)) < 0.05

    def test_weights_matter_with_three_looks(self):
        rng = np.random.default_rng(1)
        east, north, up = truth()
        asc, desc = s1_geometry(True), s1_geometry(False)
        asc2 = tuple(c * 0.98 + 0.01 for c in asc)  # a second ascending track
        clean = [los_of(g, east, north, up) for g in (asc, desc, asc2)]
        noisy = [
            clean[0] + 0.001 * rng.standard_normal((NY, NX)),
            clean[1],
            clean[2] + 0.05,
        ]
        # the biased third look is down-weighted by its large sigma
        _e_hat, u_hat, *_ = decompose_wls(
            noisy, [0.001, 0.001, 0.1], [asc, desc, asc2], north
        )
        assert np.sqrt(np.mean((u_hat - up) ** 2)) < 0.002
        _e_bad, u_bad, *_ = decompose_wls(
            noisy, [0.001, 0.001, 0.001], [asc, desc, asc2], north
        )
        assert np.sqrt(np.mean((u_bad - up) ** 2)) > 0.01

    def test_nan_looks_and_singular_geometry(self):
        east, north, up = truth()
        asc, desc = s1_geometry(True), s1_geometry(False)
        los_a = los_of(asc, east, north, up)
        los_d = los_of(desc, east, north, up)
        los_d[:, :30] = np.nan
        e_hat, _u_hat, _, _, cond = decompose_wls(
            [los_a, los_d], [0.002, 0.002], [asc, desc], north
        )
        assert np.isnan(e_hat[:, :30]).all()
        assert np.isnan(cond[:, :30]).all()
        assert np.isfinite(e_hat[:, 30:]).all()
        # two identical looks: singular
        e2, u2, _, _, _c2 = decompose_wls(
            [los_a, los_a], [0.002, 0.002], [asc, asc], north
        )
        assert np.isnan(e2).all()
        assert np.isnan(u2).all()
        # two ascending tracks: solvable but badly conditioned
        asc2 = tuple(c * 0.97 + 0.02 for c in asc)
        _, _, _, _, c3 = decompose_wls(
            [los_a, los_of(asc2, east, north, up)], [0.002, 0.002], [asc, asc2], north
        )
        assert np.nanmedian(c3) > 10

    def test_input_checks(self):
        east, north, up = truth()
        asc = s1_geometry(True)
        los_a = los_of(asc, east, north, up)
        with pytest.raises(ValueError, match="at least 2"):
            decompose_wls([los_a], [0.002], [asc], north)
        with pytest.raises(ValueError, match="one entry per look"):
            decompose_wls([los_a, los_a], [0.002], [asc, asc], north)
        with pytest.raises(ValueError, match="shape"):
            decompose_wls([los_a, los_a[:10]], [0.002, 0.002], [asc, asc], north)
        with pytest.raises(ValueError, match="negative"):
            decompose_wls([los_a, los_a], [-0.002, 0.002], [asc, asc], north)
        with pytest.raises(ValueError, match="east, north, up"):
            decompose_wls([los_a, los_a], [0.002, 0.002], [asc, asc[:2]], north)


class TestProjection:
    def test_recovers_up_with_true_horizontal_and_propagates_sigma(self):
        east, north, up = truth()
        asc = s1_geometry(True)
        los_a = los_of(asc, east, north, up)
        u_hat, su = project_vertical(
            los_a, 0.002, asc, east, north, sigma_east=0.001, sigma_north=0.0005
        )
        np.testing.assert_allclose(u_hat, up, atol=1e-12)
        e, n, u = asc
        expected = np.sqrt((0.002**2 + (e * 0.001) ** 2 + (n * 0.0005) ** 2) / u**2)
        np.testing.assert_allclose(su, expected)

    def test_unmodelled_east_leaks_into_up(self):
        """Why the mode flag exists (T53.4): short-wavelength east is absorbed as e/u."""
        east, north, up = truth()
        asc = s1_geometry(True)
        los_a = los_of(asc, east, north, up)
        e, _n, u = asc
        u_hat, _ = project_vertical(
            los_a, 0.002, asc, 0.0, north
        )  # GNSS east unknown -> 0
        leak = u_hat - up
        np.testing.assert_allclose(leak, e * east / u, atol=1e-12)
        assert np.max(np.abs(leak)) > 0.5 * np.max(np.abs(east))

    def test_low_up_pixels_are_nan(self):
        east, north, up = truth()
        asc = s1_geometry(True)
        e, n, u = asc
        u = u.copy()
        u[:, :5] = 0.05
        u_hat, su = project_vertical(
            los_of((e, n, u), east, north, up), 0.002, (e, n, u), east, north
        )
        assert np.isnan(u_hat[:, :5]).all()
        assert np.isnan(su[:, :5]).all()
        assert np.isfinite(u_hat[:, 5:]).all()


class TestDispatcher:
    def test_modes_and_fallbacks(self):
        rng = np.random.default_rng(2)
        east, north, up = truth()
        asc, desc = s1_geometry(True), s1_geometry(False)
        los_a = los_of(asc, east, north, up) + 0.002 * rng.standard_normal((NY, NX))
        los_d = los_of(desc, east, north, up) + 0.002 * rng.standard_normal((NY, NX))
        los_d[:, :40] = np.nan  # descending missing on the left
        los_a[:10, :40] = np.nan  # nothing at all in a corner
        res = decompose([los_a, los_d], [0.002, 0.002], [asc, desc], east, north)
        assert res.mode.dtype == np.uint8
        assert np.all(res.mode[10:, :40] == MODE_PROJECTION)
        assert np.all(res.mode[:10, :40] == MODE_NONE)
        assert np.all(res.mode[:, 40:] == MODE_WLS)
        assert np.isnan(res.east[:, :40]).all()
        assert np.isfinite(res.east[:, 40:]).all()
        assert np.isnan(res.up[:10, :40]).all()
        assert np.isfinite(res.up[10:, :40]).all()
        assert np.all(res.n_looks[:, 40:] == 2)
        assert np.all(res.n_looks[10:, :40] == 1)
        counts = res.counts()
        assert counts["none"] == 400
        assert counts["projection"] == 70 * 40
        # both modes are unbiased in up (projection uses the true GNSS east here)
        for sel in (res.mode == MODE_WLS, res.mode == MODE_PROJECTION):
            z = ((res.up - up) / res.sigma_up)[sel]
            assert 0.9 < np.std(z) < 1.1
            assert abs(np.mean(z)) < 0.05
        # the WLS east is consistent too
        z_e = ((res.east - east) / res.sigma_east)[res.mode == MODE_WLS]
        assert 0.9 < np.std(z_e) < 1.1

    def test_ill_conditioned_pixels_fall_back_to_projection(self):
        east, north, up = truth()
        asc = s1_geometry(True)
        asc2 = tuple(c * 0.97 + 0.02 for c in asc)
        los = [los_of(asc, east, north, up), los_of(asc2, east, north, up)]
        res = decompose(los, [0.002, 0.002], [asc, asc2], east, north, cond_max=10.0)
        assert np.all(res.mode == MODE_PROJECTION)
        assert np.nanmedian(res.cond) > 10
        np.testing.assert_allclose(
            res.up, up, atol=1e-9
        )  # projection with the true east
        relaxed = decompose(los, [0.002, 0.002], [asc, asc2], east, north, cond_max=1e6)
        assert np.all(relaxed.mode == MODE_WLS)
        # the ill-conditioned solve is right but its sigmas are large
        np.testing.assert_allclose(relaxed.up, up, atol=1e-6)
        assert np.median(relaxed.sigma_east) > 10 * np.median(res.sigma_up)

    def test_single_look_and_checks(self):
        east, north, up = truth()
        asc = s1_geometry(True)
        los_a = los_of(asc, east, north, up)
        res = decompose([los_a], [0.002], [asc], east, north)
        assert np.all(res.mode == MODE_PROJECTION)
        assert np.isnan(res.cond).all()
        np.testing.assert_allclose(res.up, up, atol=1e-12)
        with pytest.raises(ValueError, match="cond_max"):
            decompose([los_a], [0.002], [asc], east, north, cond_max=1.0)
        with pytest.raises(ValueError, match="at least 1"):
            decompose([], [], [], east, north)
