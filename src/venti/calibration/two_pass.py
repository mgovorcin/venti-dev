# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Two-pass calibration with component bookkeeping (PRD R-S1, D6; plan T33).

`calibrate_pair` turns one DISP pair plus the GNSS LOS field into a
`CalibrationResult` whose ``calibration`` is **exactly** the sum of its
components::

    calibration = cal_gnss_surface + cal_reference_offset
                + cal_tropo + cal_set + cal_unwrap_shift

so ``DISP - calibration`` is always right and every term can be audited or
undone (PRD §3.1). Two algorithms sit behind ``surface.method``:

- ``"windowed_plane"`` (gamma 0.3): delegates to `venti.surface
  .estimate_calibration_surface`, so the golden pair reproduces bit for bit;
- ``"loclin"`` (v0.5): gap filling (T29), local-linear kernel with the
  half-response cutoff (T30), coherence x robust weights (T31), remove-restore
  exclusion (T32), and optionally two passes: a robust frame-wide tie, an
  unwrap hook on ``DISP - CAL1`` whose whole-cycle shifts are applied to the
  displacement, then the final surface.

Sign convention: every component is *subtracted* from DISP. Corrections
(tropo, SET) are removed before the fit and restored into the calibration;
the unwrap shift is what was removed from DISP, so the calibrated product is
``DISP_shifted - surface``.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..spatial.resample import downsample_array, upsample_array
from ..workflow.config import CalibrationOptions
from .gaps import base_weights, fill_gaps
from .loclin import kernel_sigma_px, loclin_surface
from .remove_restore import sigma_inflation_inside
from .uncertainty import fit_sigma, resolve_k, sigma_cal
from .weights import fit_weights

logger = logging.getLogger(__name__)

__all__ = ["CalibrationResult", "UnwrapHook", "calibrate_pair"]

# ``hook(residual, valid, cycle_m) -> (shift, decisions)``: `shift` is the
# displacement to subtract (same units as the residual, 0 where no change),
# `decisions` anything the hook wants recorded (plan T36 provides it).
UnwrapHook = Callable[[np.ndarray, np.ndarray, float], tuple[np.ndarray, Any]]

COMPONENTS = (
    "cal_gnss_surface",
    "cal_reference_offset",
    "cal_tropo",
    "cal_set",
    "cal_unwrap_shift",
)


@dataclass
class CalibrationResult:
    """Everything `calibrate_pair` produced, component by component.

    All component arrays share the displacement's shape and units. A
    component that was not applied is all zeros and listed as False in
    `components_applied`.
    """

    cal_gnss_surface: np.ndarray
    cal_reference_offset: np.ndarray
    cal_tropo: np.ndarray
    cal_set: np.ndarray
    cal_unwrap_shift: np.ndarray
    components_applied: dict[str, bool]
    method: str
    reference_point: tuple[int, int]
    reference_value: float
    sigma_cal: np.ndarray | None = None
    sigma_fit: np.ndarray | None = None
    n_eff: np.ndarray | None = None
    coverage: np.ndarray | None = None
    fit_residual_std: float | None = None
    unwrap_decisions: Any = None
    n_excluded_pixels: int = 0
    passes: int = 1
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def calibration(self) -> np.ndarray:
        """The layer users subtract from DISP: the exact sum of the components."""
        total = np.zeros_like(self.cal_gnss_surface)
        for name in COMPONENTS:
            total = total + getattr(self, name)
        return total

    def assert_closed(self, atol: float = 1e-6) -> None:
        """Check ``calibration == sum(components)`` (always true by construction)."""
        total = sum(getattr(self, name) for name in COMPONENTS)
        diff = np.nanmax(np.abs(self.calibration - total))
        if not diff <= atol:
            msg = f"calibration differs from the component sum by {diff:g}"
            raise AssertionError(msg)

    def components(self) -> dict[str, np.ndarray]:
        """Return the component arrays by name."""
        return {name: getattr(self, name) for name in COMPONENTS}


def _check_shapes(disp: np.ndarray, **arrays: np.ndarray | None) -> None:
    for name, arr in arrays.items():
        if arr is not None and np.shape(arr) != disp.shape:
            msg = f"{name} shape {np.shape(arr)} != disp shape {disp.shape}"
            raise ValueError(msg)


def _reference_offset(disp: np.ndarray, ref_point: tuple[int, int]) -> float:
    refy, refx = ref_point
    value = float(disp[refy, refx])
    if not np.isfinite(value):
        msg = (
            f"Reference pixel {ref_point} has no valid displacement (value {value}); "
            "choose a reference point valid in every epoch"
        )
        raise ValueError(msg)
    return value


def _downsample_bool(mask: np.ndarray, factor: int) -> np.ndarray:
    """Downsample a mask: a block is valid when more than half of its pixels are."""
    return downsample_array(mask.astype(np.float32), factor, method="mean") > 0.5


def _upsample_nan(arr: np.ndarray, shape: tuple[int, int], factor: int) -> np.ndarray:
    """Upsample a map that may hold NaN (nearest: each pixel takes its block's node).

    `downsample_array` trims the frame to a multiple of `factor`, so the last
    rows/columns of the full grid fall beyond the last block; they take the
    nearest node rather than being cut off.
    """
    a = np.asarray(arr, dtype=np.float64)
    if factor <= 1:
        return a
    rows = np.minimum(np.arange(shape[0]) // factor, a.shape[0] - 1)
    cols = np.minimum(np.arange(shape[1]) // factor, a.shape[1] - 1)
    return a[rows][:, cols]


def _frame_cutoff_m(shape: tuple[int, int], pixel_m: float) -> float:
    """Cutoff for the pass-1 tie: four frame widths, i.e. a near-planar fit."""
    return 4.0 * max(shape) * pixel_m


# The pass-1 tie is a near-planar fit, so it is computed on a grid of at most
# this many pixels across and interpolated back: the Gaussian moments of a
# 4-frame-wide kernel on a 1300 x 1600 fit grid cost minutes, on 130 x 160
# they cost nothing, and the two agree to well below a millimetre.
TIE_MAX_PX = 256


def calibrate_pair(
    disp: np.ndarray,
    gnss_los: np.ndarray,
    mask: np.ndarray,
    ref_point: tuple[int, int],
    options: CalibrationOptions,
    pixel_m: float,
    cycle_m: float,
    *,
    gnss_los_std: np.ndarray | None = None,
    coherence: np.ndarray | None = None,
    tropo: np.ndarray | None = None,
    set_correction: np.ndarray | None = None,
    exclude_mask: np.ndarray | None = None,
    sigma_tropo: np.ndarray | float | None = None,
    sigma_ref: float = 0.0,
    unwrap_hook: UnwrapHook | None = None,
    n_jobs: int = -1,
    fit_lock: Any | None = None,
) -> CalibrationResult:
    """Calibrate one pair; return the surface and every component.

    Parameters
    ----------
    disp : np.ndarray
        LOS displacement of the pair (NaN where invalid).
    gnss_los : np.ndarray
        GNSS LOS displacement for the same pair, same grid and units
        (`venti.gnss.sampling.project_field_to_los` x the pair's span).
    mask : np.ndarray
        True = valid pixel (recommended mask AND water mask).
    ref_point : (row, col)
        Reference pixel; its displacement is removed and restored as
        ``cal_reference_offset``.
    options : CalibrationOptions
        Algorithm options (schema v2). ``surface.method`` selects the path.
    pixel_m : float
        Posting of `disp` in meters (used by the loclin cutoff).
    cycle_m : float
        One unwrapping cycle in the units of `disp` (lambda/2 of LOS
        displacement; `SensorSpec.cycle_m` when `disp` is in meters).
    gnss_los_std, coherence : np.ndarray, optional
        GNSS LOS sigma (the ``sigma_grid`` term of sigma_CAL, inflated by
        ``uncertainty.k_grid``; gamma fit weights) and temporal coherence
        (loclin weights, ``weights.coherence_power``).
    tropo, set_correction : np.ndarray, optional
        Corrections to remove before the fit; restored as ``cal_tropo`` and
        ``cal_set``.
    exclude_mask : np.ndarray, optional
        Remove-restore mask (True = excluded from the fit).
    sigma_tropo : np.ndarray or float, optional
        Uncertainty of the tropospheric correction, same units as `disp`
        (R-E1 term; 0 until a model exists).
    sigma_ref : float
        Uncertainty of the reference offset (R-E1 term), default 0.
    unwrap_hook : callable, optional
        ``hook(residual, valid, cycle_m) -> (shift, decisions)``; used only
        when ``options.unwrap_error_correction`` is True. With the loclin
        method and ``surface.two_pass`` it runs on ``DISP - CAL1``; otherwise
        on the referenced displacement minus the GNSS field.
    n_jobs, fit_lock
        Passed to the gamma windowed fit.

    """
    disp = np.array(disp, dtype=np.promote_types(np.asarray(disp).dtype, np.float32))
    _check_shapes(
        disp,
        gnss_los=gnss_los,
        mask=mask,
        gnss_los_std=gnss_los_std,
        coherence=coherence,
        tropo=tropo,
        set_correction=set_correction,
        exclude_mask=exclude_mask,
        sigma_tropo=sigma_tropo if np.ndim(sigma_tropo) > 0 else None,
    )
    zeros = np.zeros(disp.shape, dtype=disp.dtype)
    applied = {
        "cal_gnss_surface": True,
        "cal_reference_offset": True,
        "cal_tropo": tropo is not None,
        "cal_set": set_correction is not None,
        "cal_unwrap_shift": False,
    }
    cal_tropo = np.asarray(tropo, dtype=disp.dtype) if tropo is not None else zeros
    cal_set = (
        np.asarray(set_correction, dtype=disp.dtype)
        if set_correction is not None
        else zeros
    )

    # 1. corrections off, reference off: the field the fit works on
    work = disp - cal_tropo - cal_set
    ref_value = _reference_offset(work, ref_point)
    work = work - ref_value
    mask = np.asarray(mask, dtype=bool)
    valid = mask & np.isfinite(work) & np.isfinite(np.asarray(gnss_los))
    exclude = (
        np.asarray(exclude_mask, dtype=bool)
        if exclude_mask is not None
        else np.zeros(disp.shape, dtype=bool)
    )
    gnss = np.asarray(np.ma.filled(gnss_los, np.nan), dtype=np.float64)
    cal_unwrap_shift = zeros.copy()
    decisions = None

    method = options.surface.method
    if method == "windowed_plane":
        # ---- gamma path, reproduced through the gamma estimator -------------
        from ..surface import estimate_calibration_surface

        if options.unwrap_error_correction:
            # gamma shifted the referenced field in place and lost the shift;
            # do the same step here and keep it as a component
            from ..unwrap import correct_region_offset

            shifted = correct_region_offset(
                input_disp=np.where(valid, work, np.nan),
                mask=valid,
                cycle_length=cycle_m,
                min_region_area=options.unwrap.min_region_area,
            )
            shifted = np.ma.filled(shifted, np.nan)
            shift = np.where(
                np.isfinite(shifted) & np.isfinite(work), work - shifted, 0.0
            )
            cal_unwrap_shift = shift.astype(disp.dtype)
            applied["cal_unwrap_shift"] = bool(np.any(cal_unwrap_shift != 0))
            work = work - cal_unwrap_shift
        gamma_opts = options.model_copy(update={"unwrap_error_correction": False})
        window_px = max(1, round(options.window_size_meters / pixel_m))
        event_mask = ~exclude if exclude.any() else None
        gamma = estimate_calibration_surface(
            work + ref_value,  # the estimator removes and restores the offset itself
            gnss,
            valid,
            ref_point,
            window_px,
            gnss_los_std=(
                gnss_los_std if options.weight_fit_by_gnss_uncertainty else None
            ),
            event_mask=event_mask,
            options=gamma_opts,
            wavelength_m=2.0 * cycle_m,
            downsample_factor=options.downsample_factor,
            downsample_method=options.downsample_method,
            downsample_weights=coherence if options.downsample_weighted else None,
            n_jobs=n_jobs,
            fit_lock=fit_lock,
        )
        surface = np.asarray(gamma.surface, dtype=disp.dtype) - ref_value
        coverage = None
        passes = 1
        fit_std = None
        sigma = sigma_fit_map = n_eff = None  # gamma: cal-disp keeps its RBF sigma
    elif method == "loclin":
        # ---- v0.5 path ---------------------------------------------------------
        factor = max(1, int(options.downsample_factor))
        pixel_ds = pixel_m * factor

        def down(arr: np.ndarray) -> np.ndarray:
            return (
                downsample_array(arr, factor, method=options.downsample_method)
                if factor > 1
                else np.asarray(arr, float)
            )

        valid_ds = _downsample_bool(valid, factor) if factor > 1 else valid
        exclude_ds = _downsample_bool(exclude, factor) if factor > 1 else exclude
        coh_ds = (
            down(np.nan_to_num(coherence, nan=0.0)) if coherence is not None else None
        )
        gnss_ds = down(np.where(np.isfinite(gnss), gnss, 0.0))

        def fit(
            residual_ds: np.ndarray,
            cutoff_m: float,
            robust: bool,
            *,
            grids: (
                tuple[np.ndarray, np.ndarray, np.ndarray | None, float] | None
            ) = None,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            f_valid, f_exclude, f_coh, f_pixel = grids or (
                valid_ds,
                exclude_ds,
                coh_ds,
                pixel_ds,
            )
            fit_valid = f_valid & ~f_exclude & np.isfinite(residual_ds)
            if not fit_valid.any():
                msg = "no valid pixel left for the fit after masking and exclusion"
                raise ValueError(msg)
            filled, filled_mask = fill_gaps(
                np.where(fit_valid, residual_ds, np.nan), fit_valid
            )
            w = base_weights(
                fit_valid, filled_mask, options.weights.filled_pixel_weight
            )
            sigma_px = kernel_sigma_px(cutoff_m, f_pixel)
            w = fit_weights(
                filled,
                w,
                sigma_px,
                coherence=f_coh,
                coherence_power=options.weights.coherence_power,
                robust=robust,
            )
            surface_ds, cov = loclin_surface(filled, w, f_pixel, cutoff_m)
            return surface_ds, cov, w, filled

        def tie(residual_ds: np.ndarray) -> np.ndarray:
            """Pass-1 robust tie; coarse grid when the fit grid is large."""
            tie_factor = max(1, int(np.ceil(max(residual_ds.shape) / TIE_MAX_PX)))
            cutoff = _frame_cutoff_m(residual_ds.shape, pixel_ds)
            if tie_factor == 1:
                return fit(residual_ds, cutoff, robust=True)[0]
            fit_valid = valid_ds & ~exclude_ds & np.isfinite(residual_ds)
            coarse = downsample_array(
                np.where(fit_valid, residual_ds, np.nan), tie_factor, method="mean"
            )
            grids = (
                np.isfinite(coarse),
                np.zeros(coarse.shape, dtype=bool),
                downsample_array(coh_ds, tie_factor) if coh_ds is not None else None,
                pixel_ds * tie_factor,
            )
            tie_coarse = fit(coarse, cutoff, robust=True, grids=grids)[0]
            return upsample_array(tie_coarse, residual_ds.shape)

        work_ds = down(np.where(valid, work, 0.0))
        residual_ds = work_ds - gnss_ds
        passes = 1
        if options.surface.two_pass:
            passes = 2
            tie_ds = tie(residual_ds)
            if options.unwrap_error_correction and unwrap_hook is not None:
                shift_ds, decisions = unwrap_hook(
                    residual_ds - tie_ds, valid_ds & ~exclude_ds, cycle_m
                )
                shift_ds = np.asarray(shift_ds, dtype=np.float64)
                if shift_ds.shape != residual_ds.shape:
                    msg = (
                        f"unwrap hook returned shape {shift_ds.shape}, expected"
                        f" {residual_ds.shape}"
                    )
                    raise ValueError(msg)
                if np.any(shift_ds != 0):
                    shift_full = (
                        np.kron(shift_ds, np.ones((factor, factor)))[
                            : disp.shape[0], : disp.shape[1]
                        ]
                        if factor > 1
                        else shift_ds
                    )
                    cal_unwrap_shift = np.where(
                        np.isfinite(work), shift_full, 0.0
                    ).astype(disp.dtype)
                    applied["cal_unwrap_shift"] = True
                    work = work - cal_unwrap_shift
                    residual_ds = residual_ds - shift_ds
        elif options.unwrap_error_correction and unwrap_hook is not None:
            shift_ds, decisions = unwrap_hook(
                residual_ds, valid_ds & ~exclude_ds, cycle_m
            )
            shift_ds = np.asarray(shift_ds, dtype=np.float64)
            if np.any(shift_ds != 0):
                shift_full = (
                    np.kron(shift_ds, np.ones((factor, factor)))[
                        : disp.shape[0], : disp.shape[1]
                    ]
                    if factor > 1
                    else shift_ds
                )
                cal_unwrap_shift = np.where(np.isfinite(work), shift_full, 0.0).astype(
                    disp.dtype
                )
                applied["cal_unwrap_shift"] = True
                work = work - cal_unwrap_shift
                residual_ds = residual_ds - shift_ds

        surface_ds, coverage_ds, w_ds, filled_ds = fit(
            residual_ds,
            options.surface.cutoff_wavelength_meters,
            options.weights.robust,
        )
        fit_ok = (w_ds > 0) & valid_ds & ~exclude_ds
        fit_std = (
            float(
                np.sqrt(
                    np.average(
                        (residual_ds - surface_ds)[fit_ok] ** 2, weights=w_ds[fit_ok]
                    )
                )
            )
            if fit_ok.any()
            else None
        )
        surface = upsample_array(surface_ds, disp.shape) if factor > 1 else surface_ds
        coverage = (
            upsample_array(coverage_ds, disp.shape) if factor > 1 else coverage_ds
        )
        surface = np.asarray(surface, dtype=disp.dtype)

        # sigma_CAL (R-E1): fit term on the fit grid, grid term at full
        # resolution, inflation inside the interpolated areas
        sigma_px_final = kernel_sigma_px(
            options.surface.cutoff_wavelength_meters, pixel_ds
        )
        sigma_fit_ds, n_eff_ds = fit_sigma(filled_ds, w_ds, sigma_px_final, surface_ds)
        sigma_fit_map = _upsample_nan(sigma_fit_ds, disp.shape, factor)
        n_eff = _upsample_nan(n_eff_ds, disp.shape, factor)
        k = resolve_k(options.uncertainty)
        if gnss_los_std is None:
            logger.warning("no gnss_los_std given: sigma_cal carries the fit term only")
        inflation = None
        if options.uncertainty.inflate_inside_areas and exclude.any():
            inflation = sigma_inflation_inside(
                exclude,
                kernel_sigma_px(options.surface.cutoff_wavelength_meters, pixel_m),
            )
        sigma = sigma_cal(
            None if gnss_los_std is None else np.ma.filled(gnss_los_std, np.nan),
            sigma_fit_map,
            k,
            sigma_tropo=sigma_tropo,
            sigma_ref=sigma_ref,
            inflation=inflation,
        )
        sigma = np.where(valid | exclude, sigma, np.nan).astype(disp.dtype)
    else:  # pragma: no cover - Literal in the schema prevents it
        msg = f"unknown surface method {method!r}"
        raise ValueError(msg)

    result = CalibrationResult(
        cal_gnss_surface=surface,
        cal_reference_offset=np.full(disp.shape, ref_value, dtype=disp.dtype),
        cal_tropo=cal_tropo,
        cal_set=cal_set,
        cal_unwrap_shift=cal_unwrap_shift,
        components_applied=applied,
        method=method,
        reference_point=(int(ref_point[0]), int(ref_point[1])),
        reference_value=ref_value,
        sigma_cal=sigma,
        sigma_fit=sigma_fit_map,
        n_eff=n_eff,
        coverage=coverage,
        fit_residual_std=fit_std,
        unwrap_decisions=decisions,
        n_excluded_pixels=int(exclude.sum()),
        passes=passes,
    )
    result.assert_closed()
    logger.info(
        "calibrate_pair(%s, %d pass%s): components %s; fit residual std %s",
        method,
        passes,
        "es" if passes > 1 else "",
        {k: v for k, v in applied.items() if v},
        f"{fit_std:.4g}" if fit_std is not None else "n/a",
    )
    if sigma is not None:
        logger.info(
            "sigma_cal median %.4g (k = %.2f, fit term median %.4g)",
            float(np.nanmedian(sigma)),
            resolve_k(options.uncertainty),
            float(np.nanmedian(sigma_fit_map)),
        )
    return result
