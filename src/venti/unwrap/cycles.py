# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Whole-cycle unwrap-error estimation per water-mask region (PRD R-U1, D7; plan T36).

Ported from the trade-study estimator ``ref_gnss_veto`` (``unwrap_prototype/
harness.py``, 14/14 on the labelled F08882 cases) and ``two_pass/
estimator_c.py``. For one pair, on the residual ``R = DISP - CAL_1`` (the
pass-1 tie removed, so the troposphere difference across the water does not
bias the estimate), regions are visited largest first and the mainland
anchors the rest:

- **jump**: median of ``R`` on the region's coherent edge (within
  `anchor_distance_meters` of anchored land) minus the median on the anchored
  coherent land within the same distance of the region, in cycles;
- **decision**: ``round(jump) == 0`` -> no shift, the region anchors others;
  ``|jump - n| < cycle_tolerance`` -> shift by ``-n`` cycles unless the GNSS
  veto blocks it; otherwise rejected (a half-cycle offset is not an
  unwrapping error) and the region does not anchor;
- **GNSS veto**: the shift must move the region toward the GNSS field,
  ``|g - n| < |g|`` with ``g`` the region's median ``R - G`` relative to the
  mainland's. The GNSS plane cannot set the cycle count (10 km block medians
  scatter 0.4-1.2 cycles) but it tells the direction; both false shifts of
  the baseline moved away from it;
- optional **inversion-residual gate**: only regions whose median
  ``timeseries_inversion_residuals`` exceeds `residual_gate_cycles` may be
  shifted (necessary, not sufficient; off by default);
- **free offsets** (non-integer shifts) are not implemented for v0.5 and
  `UnwrapOptions.free_offsets` is rejected.

The whole module is gated: `CalibrationOptions.unwrap_error_correction` is
False until trade study TS-U1 passes, and then ``cal_unwrap_shift`` is a
component of the calibration (`venti.calibration.two_pass`).
"""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage

from ..workflow.config import UnwrapOptions

logger = logging.getLogger(__name__)

__all__ = [
    "Decisions",
    "RegionDecision",
    "apply_shifts",
    "estimate_cycles",
    "make_unwrap_hook",
    "shift_field",
]


@dataclass
class RegionDecision:
    """What was measured and decided for one region."""

    region: int
    pixels: int
    coherent: int
    jump: float = math.nan
    gnss_rel: float = math.nan
    resid_med: float = math.nan
    cycles: int = 0
    accepted: bool = False
    anchor: bool = False
    method: str = "-"
    reason: str = ""


@dataclass
class Decisions:
    """All region decisions of one pair plus what the estimator saw."""

    rows: list[RegionDecision] = field(default_factory=list)
    pair: str = ""
    mainland: int = 0
    cycle_m: float = math.nan
    pixel_m: float = math.nan

    @classmethod
    def empty(cls, pair: str = "") -> Decisions:
        """Return the decisions of a disabled correction: nothing measured."""
        return cls(pair=pair)

    @property
    def shifted(self) -> list[RegionDecision]:
        return [r for r in self.rows if r.cycles != 0]

    def by_region(self) -> dict[int, RegionDecision]:
        return {r.region: r for r in self.rows}

    def to_records(self) -> list[dict[str, Any]]:
        return [{"pair": self.pair, **asdict(r)} for r in self.rows]

    def to_csv(self, path: Path | str) -> None:
        """Write one row per region (header only when nothing was measured)."""
        fields = [
            "pair",
            *(f.name for f in RegionDecision.__dataclass_fields__.values()),
        ]
        with Path(path).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for rec in self.to_records():
                w.writerow(
                    {
                        k: f"{v:.4f}" if isinstance(v, float) else v
                        for k, v in rec.items()
                    }
                )

    def summary(self) -> dict[str, Any]:
        return {
            "pair": self.pair,
            "regions": len(self.rows),
            "measured": sum(1 for r in self.rows if math.isfinite(r.jump)),
            "shifted": [
                {"region": r.region, "cycles": r.cycles, "jump": round(r.jump, 3)}
                for r in self.shifted
            ],
            "vetoed": sum(1 for r in self.rows if r.reason.startswith("vetoed")),
            "rejected": sum(1 for r in self.rows if r.reason.startswith("rejected")),
        }


def _pixels_for_area(area_km2: float, pixel_m: float) -> int:
    return max(1, math.ceil(area_km2 * 1e6 / (pixel_m * pixel_m)))


def estimate_cycles(
    residual: np.ndarray,
    labels: np.ndarray,
    gnss_los: np.ndarray,
    coherent: np.ndarray,
    cycle_m: float,
    pixel_m: float,
    options: UnwrapOptions | None = None,
    *,
    inversion_residual: np.ndarray | None = None,
    pair: str = "",
) -> Decisions:
    """Decide the whole-cycle shift of every region (see the module docstring).

    Parameters
    ----------
    residual : np.ndarray
        ``DISP - CAL_1`` in the units of `cycle_m` (NaN where invalid).
    labels : np.ndarray
        Region labels (`venti.unwrap.regions.segment_regions`), 0 = water.
    gnss_los : np.ndarray
        GNSS LOS displacement of the pair, same grid and units.
    coherent : np.ndarray
        Pixels trusted for the medians (recommended mask AND finite residual).
    cycle_m : float
        One cycle, lambda/2 of LOS displacement (`SensorSpec.cycle_m`).
    pixel_m : float
        Posting of the grids; the option areas and distances are converted
        with it (the trade-study thresholds were 1000 / 200 pixels and 400
        pixels at 30 m).
    options : UnwrapOptions, optional
        Thresholds and switches; defaults are the trade-study values.
    inversion_residual : np.ndarray, optional
        DISP ``timeseries_inversion_residuals`` (radians); needed when
        ``options.residual_gate_cycles`` is set.
    pair : str
        Label recorded in the decisions.

    Returns
    -------
    Decisions

    """
    opts = options or UnwrapOptions()
    if opts.free_offsets:
        msg = "free offsets are not implemented (PRD R-U1: whole cycles only in v0.5)"
        raise NotImplementedError(msg)
    r_arr = np.asarray(residual, dtype=np.float64)
    lab = np.asarray(labels)
    g_arr = np.asarray(gnss_los, dtype=np.float64)
    coh = np.asarray(coherent, dtype=bool) & np.isfinite(r_arr)
    for name, a in (("labels", lab), ("gnss_los", g_arr), ("coherent", coh)):
        if a.shape != r_arr.shape:
            msg = f"{name} shape {a.shape} != residual shape {r_arr.shape}"
            raise ValueError(msg)
    if cycle_m <= 0 or pixel_m <= 0:
        msg = "cycle_m and pixel_m must be positive"
        raise ValueError(msg)
    gate = opts.residual_gate_cycles
    if gate is not None and inversion_residual is None:
        msg = "residual_gate_cycles is set but no inversion_residual was given"
        raise ValueError(msg)
    res = (
        np.asarray(inversion_residual, dtype=np.float64) / (2 * np.pi)
        if inversion_residual is not None
        else None
    )
    if res is not None and res.shape != r_arr.shape:
        msg = f"inversion_residual shape {res.shape} != residual shape {r_arr.shape}"
        raise ValueError(msg)

    near = int(opts.anchor_distance_meters / pixel_m)
    min_coh = _pixels_for_area(opts.min_coherent_area_km2, pixel_m)
    min_edge = _pixels_for_area(opts.min_edge_area_km2, pixel_m)
    half = float(cycle_m)

    ids, sizes = np.unique(lab[lab > 0], return_counts=True)
    out = Decisions(pair=pair, cycle_m=half, pixel_m=float(pixel_m))
    if ids.size == 0:
        logger.info("estimate_cycles: no regions")
        return out
    order = np.argsort(-sizes, kind="stable")
    mainland = int(ids[order[0]])
    out.mainland = mainland
    fixed = lab == mainland
    shift = np.zeros(lab.shape, dtype=np.float64)
    misfit = r_arr - g_arr
    main_coh = fixed & coh
    main_off = float(np.median(misfit[main_coh])) if main_coh.any() else 0.0

    for i in order[1:]:
        k, n = int(ids[i]), int(sizes[i])
        reg = lab == k
        ncoh = int((reg & coh).sum())
        row = RegionDecision(region=k, pixels=n, coherent=ncoh)
        out.rows.append(row)
        if ncoh < min_coh:
            row.reason = f"too small ({ncoh} < {min_coh} coherent px)"
            continue
        m = reg & coh
        row.gnss_rel = float(np.median((misfit + shift)[m]) / half) - main_off / half
        if res is not None:
            r2 = res[m]
            r2 = r2[np.isfinite(r2)]
            if r2.size:
                row.resid_med = float(np.median(r2))
        ys, xs = np.nonzero(reg)
        y0, y1 = max(int(ys.min()) - near, 0), int(ys.max()) + near + 1
        x0, x1 = max(int(xs.min()) - near, 0), int(xs.max()) + near + 1
        sub = reg[y0:y1, x0:x1]
        fx = fixed[y0:y1, x0:x1]
        cs = coh[y0:y1, x0:x1]
        if fx.any():
            d_reg = ndimage.distance_transform_edt(~sub)
            d_fix = ndimage.distance_transform_edt(~fx)
            other = fx & cs & (d_reg <= near)
            edge = sub & cs & (d_fix <= near)
            if other.sum() >= min_edge and edge.sum() >= min_edge:
                rs = (r_arr + shift)[y0:y1, x0:x1]
                row.jump = float((np.median(rs[edge]) - np.median(rs[other])) / half)
                row.method = "neighbour"
        if not math.isfinite(row.jump):
            row.reason = (
                "not measured (no anchored land within"
                f" {opts.anchor_distance_meters:g} m)"
            )
            continue
        j, g = row.jump, row.gnss_rel
        nn = round(j)
        frac = abs(j - nn)
        if nn == 0:
            row.anchor = True
            row.reason = "0 cycles"
        elif gate is not None and not row.resid_med > gate:
            row.reason = (
                f"rejected (inversion residual {row.resid_med:.2f} <= {gate:g})"
            )
        elif opts.gnss_veto and frac < opts.cycle_tolerance and abs(g - nn) >= abs(g):
            row.reason = f"vetoed by GNSS (misfit {g:+.2f} -> {g - nn:+.2f} cycles)"
        elif frac < opts.cycle_tolerance:
            row.cycles = -nn
            row.accepted = True
            row.anchor = True
            row.reason = f"shift {-nn:+d} cycles"
            shift[reg] = -nn * half
        else:
            row.reason = (
                f"rejected (|jump - {nn:d}| = {frac:.2f} >= {opts.cycle_tolerance:g})"
            )
        if row.anchor:
            fixed |= reg
    logger.info(
        "estimate_cycles%s: %d regions, %d measured, %d shifted, %d vetoed",
        f" [{pair}]" if pair else "",
        len(out.rows),
        out.summary()["measured"],
        len(out.shifted),
        out.summary()["vetoed"],
    )
    return out


def shift_field(labels: np.ndarray, decisions: Decisions, cycle_m: float) -> np.ndarray:
    """Return the displacement *added* to DISP: ``cycles * cycle_m`` per region."""
    lab = np.asarray(labels)
    out = np.zeros(lab.shape, dtype=np.float64)
    for r in decisions.shifted:
        out[lab == r.region] = r.cycles * cycle_m
    return out


def apply_shifts(
    disp: np.ndarray, labels: np.ndarray, decisions: Decisions, cycle_m: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(disp_shifted, cal_unwrap_shift)``.

    ``disp_shifted = disp + cycles * cycle_m`` on each shifted region and
    ``cal_unwrap_shift = -cycles * cycle_m`` is the component the calibration
    carries, so ``disp - cal_unwrap_shift == disp_shifted`` (D6: every term
    removed from DISP is restored into CAL). NaN pixels stay NaN and carry 0.
    """
    d = np.asarray(disp)
    if np.shape(labels) != d.shape:
        msg = f"labels shape {np.shape(labels)} != disp shape {d.shape}"
        raise ValueError(msg)
    added = shift_field(labels, decisions, cycle_m)
    added = np.where(np.isfinite(d), added, 0.0)
    shifted = d + added
    return shifted.astype(d.dtype, copy=False), (-added).astype(d.dtype, copy=False)


def make_unwrap_hook(
    labels: np.ndarray,
    gnss_los: np.ndarray,
    coherent: np.ndarray,
    options: UnwrapOptions,
    pixel_m: float,
    *,
    inversion_residual: np.ndarray | None = None,
    pair: str = "",
) -> Any:
    """Build the `venti.calibration.two_pass.UnwrapHook` for one pair.

    The arrays must be on the grid the hook will see: the fit grid of
    `calibrate_pair` (`venti.unwrap.regions.downsample_labels` and
    `venti.spatial.resample.downsample_array` with the same factor). The hook
    returns ``(shift_to_subtract, decisions)`` with
    ``shift_to_subtract = -cycles * cycle_m`` per region, which
    `calibrate_pair` records as ``cal_unwrap_shift``.
    """
    lab = np.asarray(labels)
    g_arr = np.asarray(gnss_los, dtype=np.float64)
    coh = np.asarray(coherent, dtype=bool)
    res = None if inversion_residual is None else np.asarray(inversion_residual)

    def hook(
        residual: np.ndarray, valid: np.ndarray, cycle_m: float
    ) -> tuple[np.ndarray, Decisions]:
        r = np.asarray(residual, dtype=np.float64)
        if r.shape != lab.shape:
            msg = (
                f"unwrap hook grids {lab.shape} do not match the fit grid {r.shape}; "
                "downsample labels/GNSS/coherence with the same factor"
            )
            raise ValueError(msg)
        decisions = estimate_cycles(
            r,
            lab,
            g_arr,
            coh & np.asarray(valid, dtype=bool),
            cycle_m,
            pixel_m,
            options,
            inversion_residual=res,
            pair=pair,
        )
        return -shift_field(lab, decisions, cycle_m), decisions

    return hook
