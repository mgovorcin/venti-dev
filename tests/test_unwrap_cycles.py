# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Whole-cycle estimation, shifts and the two-pass hook (PRD R-U1; plan T36.2-T36.5).

The 14-case bench on the F08882 caches (T36.4) is in
``test_unwrap_bench.py`` and needs ``VENTI_UNWRAP_BENCH_DIR``.
"""

from __future__ import annotations

import csv
import math

import numpy as np
import pytest

from venti.calibration.two_pass import calibrate_pair
from venti.unwrap.cycles import (
    Decisions,
    apply_shifts,
    estimate_cycles,
    make_unwrap_hook,
    shift_field,
)
from venti.unwrap.regions import downsample_labels, segment_regions
from venti.workflow.config import CalibrationOptions, UnwrapOptions

CYCLE = 0.02773  # m
PIXEL_M = 300.0  # anchor distance 12 km = 40 px; min coherent 0.9 km2 = 10 px


def scene():
    """Mainland + five islands with known offsets (metres).

    A: +2 cycles (real error) -> shift -2
    B: +0.6 cycle (half-cycle offset) -> rejected
    C: +1 cycle but GNSS says the island really is 1 cycle up -> vetoed
    D: +1 cycle, > 12 km from any anchored land -> not measured
    E: +1 cycle relative to A's corrected level, within 12 km of A only
       -> measured against A after A is shifted
    """
    rng = np.random.default_rng(0)
    ny, nx = 240, 360
    land = np.zeros((ny, nx), bool)
    land[:, :100] = True  # mainland
    boxes = {
        "A": (np.s_[20:60, 120:160], 2.0),
        "B": (np.s_[80:120, 120:160], 0.6),
        "C": (np.s_[140:180, 120:160], 1.0),
        "D": (np.s_[20:60, 300:350], 1.0),
        "E": (
            np.s_[20:60, 180:220],
            1.0,
        ),  # 20 px (6 km) east of A, 80 px off the mainland
    }
    for sl, _ in boxes.values():
        land[sl] = True
    labels = segment_regions(land)
    truth = np.zeros((ny, nx))
    for sl, cyc in boxes.values():
        truth[sl] = cyc * CYCLE
    residual = 0.001 * rng.standard_normal((ny, nx)) + truth
    residual[~land] = np.nan
    gnss = np.zeros((ny, nx))
    gnss[boxes["C"][0]] = 1.0 * CYCLE  # C is genuinely one cycle up per GNSS
    coherent = land.copy()
    ids = {name: int(labels[sl][20, 20]) for name, (sl, _) in boxes.items()}
    return residual, labels, gnss, coherent, ids, boxes


class TestEstimate:
    def test_decisions_on_the_scene(self):
        residual, labels, gnss, coh, ids, _ = scene()
        dec = estimate_cycles(residual, labels, gnss, coh, CYCLE, PIXEL_M, pair="t")
        by = dec.by_region()
        assert dec.mainland == int(labels[100, 50])
        a, b, c, d, e = (by[ids[k]] for k in "ABCDE")
        assert a.cycles == -2
        assert a.accepted
        assert a.anchor
        assert a.jump == pytest.approx(2.0, abs=0.05)
        assert b.cycles == 0
        assert not b.accepted
        assert not b.anchor
        assert b.reason.startswith("rejected")
        assert c.cycles == 0
        assert c.reason.startswith("vetoed by GNSS")
        assert c.gnss_rel == pytest.approx(0.0, abs=0.1)
        assert d.cycles == 0
        assert not math.isfinite(d.jump)
        assert d.reason.startswith("not measured")
        # E is measured against A's shifted level
        assert e.cycles == -1
        assert e.method == "neighbour"
        assert e.jump == pytest.approx(1.0, abs=0.05)
        assert [r.region for r in dec.shifted] == sorted(
            [ids["A"], ids["E"]], key=lambda k: -by[k].pixels
        )
        s = dec.summary()
        assert s["measured"] == 4
        assert len(s["shifted"]) == 2
        assert s["vetoed"] == 1

    def test_veto_off_shifts_c_and_gate_needs_residual(self):
        residual, labels, gnss, coh, ids, _ = scene()
        dec = estimate_cycles(
            residual, labels, gnss, coh, CYCLE, PIXEL_M, UnwrapOptions(gnss_veto=False)
        )
        assert dec.by_region()[ids["C"]].cycles == -1
        gated = UnwrapOptions(residual_gate_cycles=0.05)
        with pytest.raises(ValueError, match="inversion_residual"):
            estimate_cycles(residual, labels, gnss, coh, CYCLE, PIXEL_M, gated)
        res = np.zeros(residual.shape)
        res[labels == ids["A"]] = 2 * np.pi * 0.4  # 0.4 cycles misfit on A only
        dec = estimate_cycles(
            residual, labels, gnss, coh, CYCLE, PIXEL_M, gated, inversion_residual=res
        )
        by = dec.by_region()
        assert by[ids["A"]].cycles == -2
        assert by[ids["A"]].resid_med == pytest.approx(0.4)
        assert by[ids["E"]].cycles == 0
        assert "inversion residual" in by[ids["E"]].reason

    def test_zero_cycle_regions_anchor_and_small_ones_do_not(self):
        residual, labels, gnss, coh, ids, boxes = scene()
        residual[boxes["A"][0]] -= 2 * CYCLE  # A is now at 0 (+ noise)
        residual[boxes["A"][0]] += 0.35 * CYCLE  # ... but 0.35 cycles off
        residual[boxes["E"][0]] += 0.35 * CYCLE  # E keeps its +1 cycle against A
        dec = estimate_cycles(residual, labels, gnss, coh, CYCLE, PIXEL_M)
        by = dec.by_region()
        assert by[ids["A"]].cycles == 0
        assert by[ids["A"]].anchor
        assert by[ids["E"]].cycles == -1  # E still measured through A
        coh[boxes["A"][0]] = False
        dec = estimate_cycles(residual, labels, gnss, coh, CYCLE, PIXEL_M)
        by = dec.by_region()
        assert by[ids["A"]].reason.startswith("too small")
        assert not by[ids["A"]].anchor
        assert by[ids["E"]].reason.startswith("not measured")

    def test_inputs_are_checked(self):
        residual, labels, gnss, coh, *_ = scene()
        with pytest.raises(ValueError, match="labels shape"):
            estimate_cycles(residual, labels[:10], gnss, coh, CYCLE, PIXEL_M)
        with pytest.raises(ValueError, match="positive"):
            estimate_cycles(residual, labels, gnss, coh, 0.0, PIXEL_M)
        with pytest.raises(NotImplementedError, match="free offsets"):
            estimate_cycles(
                residual,
                labels,
                gnss,
                coh,
                CYCLE,
                PIXEL_M,
                UnwrapOptions(free_offsets=True),
            )
        empty = estimate_cycles(
            residual, np.zeros_like(labels), gnss, coh, CYCLE, PIXEL_M
        )
        assert empty.rows == []
        assert empty.mainland == 0


class TestShifts:
    def test_apply_shifts_closes_with_the_component(self):
        residual, labels, gnss, coh, ids, boxes = scene()
        dec = estimate_cycles(residual, labels, gnss, coh, CYCLE, PIXEL_M)
        disp = residual.astype(np.float32)
        shifted, comp = apply_shifts(disp, labels, dec, CYCLE)
        assert shifted.dtype == comp.dtype == np.float32
        np.testing.assert_allclose(
            np.nan_to_num(disp - comp), np.nan_to_num(shifted), atol=1e-8
        )
        a = boxes["A"][0]
        np.testing.assert_allclose(comp[a], 2 * CYCLE, rtol=1e-5)  # subtracted
        assert np.nanmedian(shifted[a]) == pytest.approx(0.0, abs=0.001)
        assert np.all(comp[labels == ids["B"]] == 0)
        assert np.isnan(shifted[~np.isfinite(disp)]).all()
        assert np.all(comp[~np.isfinite(disp)] == 0)
        np.testing.assert_allclose(shift_field(labels, dec, CYCLE), -comp, atol=1e-8)
        with pytest.raises(ValueError, match="labels shape"):
            apply_shifts(disp, labels[:3], dec, CYCLE)

    def test_decisions_csv_round_trip_and_empty(self, tmp_path):
        residual, labels, gnss, coh, ids, _ = scene()
        dec = estimate_cycles(residual, labels, gnss, coh, CYCLE, PIXEL_M, pair="p1")
        path = tmp_path / "decisions.csv"
        dec.to_csv(path)
        with path.open() as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == len(dec.rows) == 5
        a = next(r for r in rows if int(r["region"]) == ids["A"])
        assert a["pair"] == "p1"
        assert int(a["cycles"]) == -2
        assert a["accepted"] == "True"
        empty = Decisions.empty("p2")
        empty.to_csv(path)
        with path.open() as fh:
            text = fh.read().splitlines()
        assert len(text) == 1
        assert text[0].startswith("pair,region,")
        assert empty.shifted == []
        assert empty.summary()["shifted"] == []


class TestHook:
    """T36.5: the estimator through `calibrate_pair`'s hook."""

    @staticmethod
    def _scene():
        """Two-pass test scene: plane + 1-cycle island behind a channel."""
        rng = np.random.default_rng(3)
        ny, nx = 150, 180
        pixel_m = 180.0
        yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
        gnss = 0.004 * xx / nx - 0.002 * yy / ny
        ramp = 0.010 * (xx + yy) / (nx + ny) + 0.003
        land = np.ones((ny, nx), bool)
        land[:, :8] = False
        land[100:150, 120:180] = False
        island = np.zeros((ny, nx), bool)
        island[110:145, 135:175] = True
        land |= island
        disp = gnss + ramp + 0.0005 * rng.standard_normal((ny, nx))
        disp[island] += CYCLE
        disp[~land] = np.nan
        return disp, gnss, land, island, pixel_m

    def test_hook_sees_disp_minus_cal1_with_the_gnss_signal(self):
        """Regression: the hook gets DISP - CAL1, not DISP - GNSS - CAL1.

        `estimate_cycles` subtracts the GNSS field itself for the veto; given
        a residual with GNSS already removed it subtracted it twice. With
        DISP equal to the GNSS field the pass-1 tie is ~0, so the hook must
        see the GNSS field (before the fix it saw ~0).
        """
        ny, nx = 90, 120
        yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
        gnss = 0.02 * xx / nx - 0.01 * yy / ny
        disp = gnss.copy()
        land = np.ones((ny, nx), bool)
        seen = {}

        def spy(residual, valid, cycle_m):
            seen["residual"], seen["valid"] = residual, valid
            return np.zeros_like(residual), None

        opts = CalibrationOptions(
            unwrap_error_correction=True,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 12_000.0,
                "two_pass": True,
            },
        )
        calibrate_pair(disp, gnss, land, (0, 0), opts, 180.0, CYCLE, unwrap_hook=spy)
        r = seen["residual"][seen["valid"]]
        expected = gnss - gnss[0, 0]  # the reference pixel is removed
        np.testing.assert_allclose(
            seen["residual"][seen["valid"]], expected[seen["valid"]], atol=2e-4
        )
        assert np.ptp(r) > 0.02  # the GNSS signal is there

    def test_hook_removes_the_island_cycle(self):
        disp, gnss, land, island, pixel_m = self._scene()
        labels = segment_regions(land)
        opts = CalibrationOptions(
            unwrap_error_correction=True,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 12_000.0,
                "two_pass": True,
            },
            unwrap={"min_coherent_area_km2": 0.3, "min_edge_area_km2": 0.1},
        )
        hook = make_unwrap_hook(labels, gnss, land, opts.unwrap, pixel_m, pair="x")
        res = calibrate_pair(
            disp, gnss, land, (20, 100), opts, pixel_m, CYCLE, unwrap_hook=hook
        )
        assert res.components_applied["cal_unwrap_shift"]
        np.testing.assert_allclose(res.cal_unwrap_shift[island], CYCLE, rtol=1e-5)
        assert np.all(res.cal_unwrap_shift[~island] == 0)
        assert isinstance(res.unwrap_decisions, Decisions)
        assert res.unwrap_decisions.pair == "x"
        assert [r.cycles for r in res.unwrap_decisions.shifted] == [-1]
        calibrated = disp - res.calibration
        err = (calibrated - gnss)[land]
        assert np.sqrt(np.mean(err**2)) < 0.001
        res.assert_closed()

    def test_hook_on_the_downsampled_grid(self):
        disp, gnss, land, _island, pixel_m = self._scene()
        labels = segment_regions(land)
        factor = 3
        from venti.spatial.resample import downsample_array

        small_labels = downsample_labels(labels, factor)
        small_gnss = downsample_array(gnss, factor)
        small_land = downsample_array(land.astype(np.float32), factor) > 0.5
        opts = CalibrationOptions(
            unwrap_error_correction=True,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 12_000.0,
                "two_pass": True,
            },
            unwrap={"min_coherent_area_km2": 0.3, "min_edge_area_km2": 0.1},
        )
        opts.downsample_factor = factor
        hook = make_unwrap_hook(
            small_labels, small_gnss, small_land, opts.unwrap, pixel_m * factor
        )
        res = calibrate_pair(
            disp, gnss, land, (20, 100), opts, pixel_m, CYCLE, unwrap_hook=hook
        )
        assert res.components_applied["cal_unwrap_shift"]
        # interior of the island: one cycle; block edges may be mixed
        assert np.median(res.cal_unwrap_shift[115:140, 140:170]) == pytest.approx(
            CYCLE, rel=1e-5
        )
        # a hook on the wrong grid is refused
        bad = make_unwrap_hook(labels, gnss, land, opts.unwrap, pixel_m)
        with pytest.raises(ValueError, match="do not match the fit grid"):
            calibrate_pair(
                disp, gnss, land, (20, 100), opts, pixel_m, CYCLE, unwrap_hook=bad
            )

    def test_hook_for_pair_builds_the_fit_grid_itself(self):
        """`unwrap_hook_for_pair` from full-resolution layers, factor 3."""
        from venti.unwrap import unwrap_hook_for_pair

        disp, gnss, land, _island, pixel_m = self._scene()
        opts = CalibrationOptions(
            unwrap_error_correction=True,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 12_000.0,
                "two_pass": True,
            },
            unwrap={"min_coherent_area_km2": 0.3, "min_edge_area_km2": 0.1},
        )
        opts.downsample_factor = 3
        hook = unwrap_hook_for_pair(
            land,
            np.isfinite(disp),
            gnss,
            land & np.isfinite(disp),
            opts,
            pixel_m,
            pair="p",
        )
        res = calibrate_pair(
            disp, gnss, land, (20, 100), opts, pixel_m, CYCLE, unwrap_hook=hook
        )
        assert [r.cycles for r in res.unwrap_decisions.shifted] == [-1]
        assert res.unwrap_decisions.pair == "p"
        assert np.median(res.cal_unwrap_shift[115:140, 140:170]) == pytest.approx(
            CYCLE, rel=1e-5
        )
        assert np.all(res.cal_unwrap_shift[:100] == 0)
        res.assert_closed()

    def test_shift_on_a_frame_not_divisible_by_the_factor(self):
        """Regression: a non-zero shift on a frame whose size is not a multiple
        of the downsample factor raised a broadcast error (every shifted epoch
        of the F08882 e2e failed: 7733 x 9464 vs a 7728 x 9462 shift)."""
        from venti.unwrap import unwrap_hook_for_pair

        disp, gnss, land, _island, pixel_m = self._scene()
        # 148 x 176: blocks of 3 cover 147 x 174, the island reaches column 174
        disp, gnss, land = disp[:-2, :-4], gnss[:-2, :-4], land[:-2, :-4]
        opts = CalibrationOptions(
            unwrap_error_correction=True,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 12_000.0,
                "two_pass": True,
            },
            unwrap={"min_coherent_area_km2": 0.3, "min_edge_area_km2": 0.1},
        )
        opts.downsample_factor = 3
        hook = unwrap_hook_for_pair(
            land, np.isfinite(disp), gnss, land & np.isfinite(disp), opts, pixel_m
        )
        res = calibrate_pair(
            disp, gnss, land, (20, 100), opts, pixel_m, CYCLE, unwrap_hook=hook
        )
        assert res.cal_unwrap_shift.shape == disp.shape
        assert [r.cycles for r in res.unwrap_decisions.shifted] == [-1]
        # the trailing rows/columns beyond the last block are covered too
        assert res.cal_unwrap_shift[130, 174] == pytest.approx(CYCLE, rel=1e-5)
        res.assert_closed()

    def test_disabled_means_zero_component_and_no_decisions(self):
        disp, gnss, land, _island, pixel_m = self._scene()
        labels = segment_regions(land)
        opts = CalibrationOptions(
            unwrap_error_correction=False,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 12_000.0,
                "two_pass": True,
            },
        )
        hook = make_unwrap_hook(labels, gnss, land, opts.unwrap, pixel_m)
        res = calibrate_pair(
            disp, gnss, land, (20, 100), opts, pixel_m, CYCLE, unwrap_hook=hook
        )
        assert not res.components_applied["cal_unwrap_shift"]
        assert np.all(res.cal_unwrap_shift == 0)
        assert res.unwrap_decisions is None
