# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tropospheric correction modes (PRD R-T1; plan T34).

Two layers:

1. **The cal-disp numerics, ported verbatim.** `interpolate_in_time`,
   `interpolate_to_dem_surface` (zenith delay cube -> DEM surface, 512-row
   blocks) and `compute_los_correction` (``-ZTD / los_up``) are copies of
   ``cal_disp.product._tropo`` without the file I/O; `pair_correction` is the
   ``sec - ref`` the cal-disp workflow forms. ``tests/test_calibration_tropo.py``
   keeps them bit-identical to cal-disp while cal-disp still carries its own
   copy (plan T37 switches it to these).
2. **Modes.** Trade study ``tropo_check`` (four frames, 2026-10) found that
   the full HRES correction costs calibration-residual sill on every frame up
   to 0.4 km relief and that only its height-dependent part has skill where
   there is relief (Los Angeles, 1.55 km: -15 % sill). Hence:

   - ``off``: no tropo component;
   - ``full``: the whole LOS delay difference (gamma behaviour);
   - ``stratified``: per epoch, the least-squares fit
     ``delay ~ a + b x + c y + d h + e h**2`` to the delay (``h`` DEM height,
     ``x``/``y`` map coordinates). The plane is absorbed by the calibration
     surface anyway; the lateral 10-50 km structure, where the model showed
     no skill, is dropped and cannot inject turbulent model noise;
   - ``auto``: from the DEM relief (p5-p95): off below
     `TropoOptions.relief_off_meters`, stratified at or above
     `TropoOptions.relief_stratified_meters`, off in between until TS-T1;
   - ``legacy``: defer to ``apply_tropo_correction`` (schema v1 behaviour).

   The frame-parameter table (plan T27) pins the mode per frame; `auto` is
   the rule for frames not in the table.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime
from typing import Literal

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from rasterio.enums import Resampling
from scipy.interpolate import RegularGridInterpolator

from ..workflow.config import TropoOptions

logger = logging.getLogger(__name__)

__all__ = [
    "StratifiedModel",
    "TropoMode",
    "apply_tropo",
    "choose_tropo_mode",
    "compute_los_correction",
    "dem_relief",
    "interpolate_in_time",
    "interpolate_to_dem_surface",
    "pair_correction",
    "stratified_delay",
    "stratified_tropo",
]

TropoMode = Literal["off", "stratified", "full"]
"""A resolved mode: what `apply_tropo` does. ``auto`` and ``legacy`` resolve to one."""

# The stratified fit subsamples the delay grid for the normal equations; the
# prototype used every 7th pixel (~4 km on a 600 m grid).  The model has five
# coefficients, so any step that leaves a few thousand pixels is equivalent.
DEFAULT_SUBSAMPLE = 7


# --------------------------------------------------------------------------
# 1. cal-disp numerics (keep identical to cal_disp.product._tropo)


def interpolate_in_time(
    da_early: xr.DataArray,
    da_late: xr.DataArray,
    early_date: datetime,
    late_date: datetime,
    target_datetime: datetime,
    *,
    early_name: str = "",
    late_name: str = "",
) -> xr.DataArray:
    """Linearly interpolate two zenith-delay cubes in time.

    Parameters
    ----------
    da_early, da_late : xr.DataArray
        Total zenith delay cubes (``height, y, x``) at `early_date` and
        `late_date` (what ``TropoProduct.get_total_delay`` returns).
    early_date, late_date : datetime
        Model times of the two cubes; `early_date` must precede `late_date`.
    target_datetime : datetime
        Acquisition time; must lie in ``[early_date, late_date]``.
    early_name, late_name : str, optional
        Product names recorded in the attributes.

    Returns
    -------
    xr.DataArray
        ``(1 - w) * early + w * late`` with ``w`` the time fraction.

    Raises
    ------
    ValueError
        If the dates are out of order or the target is outside them.

    """
    if early_date > late_date:
        raise ValueError(
            f"Early product date ({early_date}) must be before "
            f"late product date ({late_date})"
        )
    if target_datetime < early_date or target_datetime > late_date:
        raise ValueError(
            f"Target datetime ({target_datetime}) must be between "
            f"early ({early_date}) and late ({late_date}) dates"
        )

    delta_total = (late_date - early_date).total_seconds()
    delta_target = (target_datetime - early_date).total_seconds()
    weight = delta_target / delta_total

    da_interp = (1 - weight) * da_early + weight * da_late
    da_interp.name = "zenith_total_delay"
    da_interp.attrs.update(
        {
            "interpolation_method": "linear",
            "early_product": early_name,
            "late_product": late_name,
            "early_date": early_date.isoformat(),
            "late_date": late_date.isoformat(),
            "target_date": target_datetime.isoformat(),
            "interpolation_weight": float(weight),
            "long_name": "Total zenith tropospheric delay",
            "units": "meters",
        }
    )
    return da_interp


def interpolate_to_dem_surface(
    da_tropo_cube: xr.DataArray,
    dem: xr.DataArray,
    method: str = "linear",
    block_rows: int = 512,
) -> xr.DataArray:
    """Interpolate a 3D zenith-delay cube to the DEM surface heights.

    Parameters
    ----------
    da_tropo_cube : xr.DataArray
        Delay cube with dims ``(height, y, x)`` (or ``latitude``/``longitude``).
        Assumed EPSG:4326 when it carries no CRS.
    dem : xr.DataArray
        Surface heights in metres on the target grid, with a CRS
        (``dem.rio.crs``).
    method : str, optional
        ``"linear"`` or ``"nearest"``. Default ``"linear"``.
    block_rows : int, optional
        DEM rows interpolated per call; bounds the interpolator's memory and
        does not change the result. Default 512.

    Returns
    -------
    xr.DataArray
        float32 delay at the DEM surface on the DEM grid (NaN outside the cube).

    Raises
    ------
    ValueError
        If the DEM has no CRS or the cube has no ``height`` dimension.

    """
    if "latitude" in da_tropo_cube.dims:
        da_tropo_cube = da_tropo_cube.rename({"latitude": "y", "longitude": "x"})

    if not hasattr(dem, "rio") or dem.rio.crs is None:
        raise ValueError(
            "DEM is missing CRS information. "
            "Use dem.rio.write_crs() to set the CRS before calling this function."
        )
    dem_crs = dem.rio.crs

    if not hasattr(da_tropo_cube, "rio") or da_tropo_cube.rio.crs is None:
        da_tropo_cube = da_tropo_cube.rio.write_crs("EPSG:4326")

    if da_tropo_cube.rio.crs != dem_crs:
        td_utm = da_tropo_cube.rio.reproject(dem_crs, resampling=Resampling.cubic)
    else:
        td_utm = da_tropo_cube

    if "height" not in td_utm.dims:
        raise ValueError(
            f"No height dimension found. Available dims: {list(td_utm.dims)}"
        )

    rgi = RegularGridInterpolator(
        (td_utm["height"].values, td_utm.y.values, td_utm.x.values),
        np.nan_to_num(td_utm.values),
        method=method,
        bounds_error=False,
        fill_value=np.nan,
    )

    # A new float32 array on the DEM grid (coordinates and CRS kept).  Never
    # write into the DEM's own buffer: the DISP-S1-STATIC DEM is float16 on
    # disk, and a delay of ~2.4 m stored in float16 is quantised to ~2 mm.
    out = dem.astype(np.float32, copy=True)

    # Row blocks: on a full frame (~73 M pixels) a single call needs ~12 GB of
    # float64/int64 temporaries inside the interpolator.  Each point is
    # interpolated independently, so the result does not depend on the block.
    y = dem.y.values
    x = dem.x.values
    dem_values = dem.values
    out_values = out.values
    for start in range(0, dem.shape[0], block_rows):
        stop = min(start + block_rows, dem.shape[0])
        yy, xx = np.meshgrid(y[start:stop], x, indexing="ij")
        pts = np.column_stack([dem_values[start:stop].ravel(), yy.ravel(), xx.ravel()])
        out_values[start:stop] = rgi(pts).reshape(stop - start, -1).astype(np.float32)
    out.name = da_tropo_cube.name or "tropospheric_delay"
    out.attrs.update(
        {
            "interpolation_method": method,
            "interpolated_from": "3D tropospheric model",
            "units": "meters",
            "long_name": "Tropospheric delay at DEM surface",
        }
    )
    if "time" in td_utm.coords:
        out.attrs["time"] = str(td_utm.time.values)
    if "target_date" in da_tropo_cube.attrs:
        out.attrs["target_date"] = da_tropo_cube.attrs["target_date"]
    if "interpolation_weight" in da_tropo_cube.attrs:
        out.attrs["interpolation_weight"] = da_tropo_cube.attrs["interpolation_weight"]
    return out


def compute_los_correction(
    zenith_delay_2d: xr.DataArray,
    los_up: xr.DataArray,
    reference_correction: xr.DataArray | None = None,
) -> xr.DataArray:
    """Convert a zenith delay to a line-of-sight correction, ``-ZTD / los_up``.

    The sign follows the DISP convention (positive = apparent uplift towards
    the satellite). When `reference_correction` is given it is subtracted and
    the result is float32, as in cal-disp.
    """
    if los_up.shape != zenith_delay_2d.shape:
        if hasattr(los_up, "rio") and hasattr(zenith_delay_2d, "rio"):
            los_up = los_up.rio.reproject_match(zenith_delay_2d)
        else:
            raise ValueError(
                f"Shape mismatch: los_up {los_up.shape} vs "
                f"zenith_delay {zenith_delay_2d.shape}"
            )

    los_correction = -1 * (zenith_delay_2d / los_up)

    if reference_correction is not None:
        if reference_correction.shape != los_correction.shape and hasattr(
            reference_correction, "rio"
        ):
            reference_correction = reference_correction.rio.reproject_match(
                los_correction
            )
        los_correction = (los_correction - reference_correction).astype(np.float32)

    los_correction.name = "los_correction"
    los_correction.attrs.update(
        {
            "units": "meters",
            "long_name": "Line-of-sight atmospheric correction",
            "line_of_sight_convention": (
                "Positive means decrease in delay (apparent uplift towards satellite)"
            ),
        }
    )
    if reference_correction is not None:
        los_correction.attrs["reference_subtracted"] = "yes"
    return los_correction


def pair_correction(ref_los: np.ndarray, sec_los: np.ndarray) -> np.ndarray:
    """Return the pair's LOS tropo correction ``sec - ref`` as float32.

    This is what the cal-disp workflow removes before the fit (and what the
    ``full`` mode restores as ``cal_tropo``).
    """
    ref = np.asarray(ref_los, dtype=np.float32)
    sec = np.asarray(sec_los, dtype=np.float32)
    if ref.shape != sec.shape:
        msg = f"reference {ref.shape} and secondary {sec.shape} layers differ"
        raise ValueError(msg)
    return sec - ref


# --------------------------------------------------------------------------
# 2. Modes


def dem_relief(
    dem: np.ndarray,
    valid: np.ndarray | None = None,
    percentiles: tuple[float, float] = (5.0, 95.0),
) -> float:
    """Return the DEM relief in metres: the p95 - p5 height spread.

    The spread, not the range, so a handful of towers or voids does not
    decide the mode. NaN heights (and pixels where `valid` is False) are
    ignored.
    """
    h = np.asarray(dem, dtype=np.float64)
    ok = np.isfinite(h)
    if valid is not None:
        ok &= np.asarray(valid, dtype=bool)
    if ok.sum() < 2:
        msg = "dem_relief needs at least two finite heights"
        raise ValueError(msg)
    lo, hi = np.percentile(h[ok], percentiles)
    return float(hi - lo)


def _design(x: np.ndarray, y: np.ndarray, h: np.ndarray) -> np.ndarray:
    return np.column_stack([np.ones(h.size), x, y, h, h * h])


@dataclass(frozen=True)
class StratifiedModel:
    """``delay = a + b x + c y + d h + e h**2`` fitted to one delay field.

    `x`, `y` are map coordinates divided by `xy_scale` and `h` heights divided
    by `h_scale` (defaults 1e5 m and 1e3 m, so ``d`` is in delay units per km).
    """

    coefficients: np.ndarray
    xy_scale: float = 1e5
    h_scale: float = 1e3
    n_fit: int = 0

    @property
    def gradient_per_km(self) -> float:
        """Linear height coefficient ``d`` in delay units per km of height."""
        return float(self.coefficients[3]) * 1e3 / self.h_scale

    def evaluate(
        self, dem: np.ndarray, x: np.ndarray, y: np.ndarray, valid: np.ndarray
    ) -> np.ndarray:
        """Return the model on the grid (NaN where `valid` is False)."""
        out = np.full(dem.shape, np.nan)
        h = np.asarray(dem, dtype=np.float64)
        out[valid] = (
            _design(
                x[valid] / self.xy_scale, y[valid] / self.xy_scale, h[valid] / self.h_scale
            )
            @ self.coefficients
        )
        return out


def _grid_xy(
    shape: tuple[int, int], x: np.ndarray | None, y: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = shape
    if x is None or y is None:
        if (x is None) != (y is None):
            msg = "give both x and y or neither"
            raise ValueError(msg)
        yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
        return xx, yy
    xa = np.asarray(x, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64)
    if xa.ndim == 1 and ya.ndim == 1:
        if xa.size != nx or ya.size != ny:
            msg = f"x ({xa.size}) and y ({ya.size}) do not match the grid {shape}"
            raise ValueError(msg)
        return np.meshgrid(xa, ya)
    if xa.shape != shape or ya.shape != shape:
        msg = f"x {xa.shape} / y {ya.shape} do not match the grid {shape}"
        raise ValueError(msg)
    return xa, ya


def stratified_tropo(
    delay: np.ndarray,
    dem: np.ndarray,
    valid: np.ndarray | None = None,
    *,
    x: np.ndarray | None = None,
    y: np.ndarray | None = None,
    subsample: int = DEFAULT_SUBSAMPLE,
    xy_scale: float = 1e5,
    h_scale: float = 1e3,
) -> StratifiedModel:
    """Fit the stratified (height-dependent) model to one LOS delay field.

    Parameters
    ----------
    delay : np.ndarray
        LOS tropospheric delay (or a pair difference), any units.
    dem : np.ndarray
        Heights in metres on the same grid.
    valid : np.ndarray, optional
        Pixels to use (True = use); non-finite pixels are dropped anyway.
    x, y : np.ndarray, optional
        Map coordinates (1-D axes or 2-D grids) in metres; pixel indices
        when omitted, with `xy_scale` then meaningless but harmless (the
        fitted field does not depend on the parametrisation).
    subsample : int
        Use every `subsample`-th row and column for the normal equations.
    xy_scale, h_scale : float
        Conditioning scales for the design matrix.

    Returns
    -------
    StratifiedModel

    Raises
    ------
    ValueError
        If shapes differ, or fewer than 50 pixels remain for five unknowns.

    """
    d = np.asarray(delay, dtype=np.float64)
    h = np.asarray(dem, dtype=np.float64)
    if h.shape != d.shape:
        msg = f"DEM {h.shape} not on the delay grid {d.shape}"
        raise ValueError(msg)
    if subsample < 1:
        msg = "subsample must be >= 1"
        raise ValueError(msg)
    ok = np.isfinite(d) & np.isfinite(h)
    if valid is not None:
        ok &= np.asarray(valid, dtype=bool)
    xx, yy = _grid_xy(d.shape, x, y)
    s = (slice(None, None, subsample), slice(None, None, subsample))
    m = ok[s]
    n_fit = int(m.sum())
    if n_fit < 50:
        msg = f"stratified fit needs >= 50 pixels after subsampling, got {n_fit}"
        raise ValueError(msg)
    a = _design(xx[s][m] / xy_scale, yy[s][m] / xy_scale, h[s][m] / h_scale)
    coef = np.linalg.lstsq(a, d[s][m], rcond=None)[0]
    model = StratifiedModel(coef, xy_scale, h_scale, n_fit)
    logger.info(
        "stratified tropo: d=%+.2f mm/km, e=%+.2f mm/km^2 from %d pixels",
        1e3 * model.gradient_per_km,
        1e3 * float(coef[4]) * (1e3 / h_scale) ** 2,
        n_fit,
    )
    return model


def stratified_delay(
    delay: np.ndarray,
    dem: np.ndarray,
    valid: np.ndarray | None = None,
    *,
    x: np.ndarray | None = None,
    y: np.ndarray | None = None,
    subsample: int = DEFAULT_SUBSAMPLE,
) -> np.ndarray:
    """Return the stratified part of `delay` on the full grid (NaN where invalid).

    Linear in `delay`, so the stratified pair correction equals the difference
    of the per-date stratified layers whenever both dates share the valid set.
    """
    d = np.asarray(delay, dtype=np.float64)
    h = np.asarray(dem, dtype=np.float64)
    model = stratified_tropo(d, h, valid, x=x, y=y, subsample=subsample)
    ok = np.isfinite(d) & np.isfinite(h)
    if valid is not None:
        ok &= np.asarray(valid, dtype=bool)
    xx, yy = _grid_xy(d.shape, x, y)
    return model.evaluate(h, xx, yy, ok)


def choose_tropo_mode(
    options: TropoOptions,
    dem: np.ndarray | None = None,
    *,
    apply_tropo_correction: bool = True,
    valid: np.ndarray | None = None,
) -> TropoMode:
    """Resolve `TropoOptions.mode` to what `apply_tropo` does.

    Parameters
    ----------
    options : TropoOptions
        The frame's tropo options, i.e. ``calibration_options.tropo`` after
        the frame-parameter table (plan T27) has been applied.
    dem : np.ndarray, optional
        Heights in metres; required for ``auto``.
    apply_tropo_correction : bool
        The schema-v1 switch, consulted for ``legacy`` only.
    valid : np.ndarray, optional
        Pixels counted in the relief (e.g. the frame's land mask).

    Returns
    -------
    {"off", "stratified", "full"}

    """
    mode = options.mode
    if mode == "legacy":
        return "full" if apply_tropo_correction else "off"
    if mode in ("off", "stratified", "full"):
        return mode  # type: ignore[return-value]
    if dem is None:
        msg = "tropo mode 'auto' needs the DEM to measure the relief"
        raise ValueError(msg)
    relief = dem_relief(dem, valid)
    if relief >= options.relief_stratified_meters:
        chosen: TropoMode = "stratified"
    else:
        # Below the stratified threshold the correction costs sill (tropo_check:
        # -7 % to -10 % on < 0.3 km frames, +5-8 % sill at 0.41 km); the band
        # up to relief_stratified_meters stays off until TS-T1 says otherwise.
        chosen = "off"
    logger.info(
        "tropo mode auto: relief p5-p95 = %.0f m -> %s (off < %.0f m, "
        "stratified >= %.0f m)",
        relief,
        chosen,
        options.relief_off_meters,
        options.relief_stratified_meters,
    )
    return chosen


def apply_tropo(
    correction: np.ndarray | None,
    mode: TropoMode,
    *,
    dem: np.ndarray | None = None,
    valid: np.ndarray | None = None,
    x: np.ndarray | None = None,
    y: np.ndarray | None = None,
    subsample: int = DEFAULT_SUBSAMPLE,
) -> np.ndarray | None:
    """Return the ``cal_tropo`` component for a resolved mode.

    Parameters
    ----------
    correction : np.ndarray or None
        The pair's full LOS correction (`pair_correction`), or None when no
        tropo products were staged.
    mode : {"off", "stratified", "full"}
        A resolved mode (`choose_tropo_mode`).
    dem, valid, x, y, subsample
        Passed to `stratified_delay` for ``stratified``.

    Returns
    -------
    np.ndarray or None
        None for ``off`` (component not applied: `calibrate_pair` records
        ``components_applied["cal_tropo"] = False``), the full field for
        ``full``, its stratified part for ``stratified``. The dtype of
        `correction` is kept.

    Raises
    ------
    ValueError
        For ``full``/``stratified`` without a correction, or ``stratified``
        without a DEM.

    """
    if mode == "off":
        if correction is not None:
            logger.info("tropo mode off: ignoring the staged tropo correction")
        return None
    if correction is None:
        msg = f"tropo mode {mode!r} needs the tropo correction, got None"
        raise ValueError(msg)
    corr = np.asarray(correction)
    if mode == "full":
        return corr
    if mode == "stratified":
        if dem is None:
            msg = "tropo mode 'stratified' needs the DEM"
            raise ValueError(msg)
        strat = stratified_delay(corr, dem, valid, x=x, y=y, subsample=subsample)
        return strat.astype(corr.dtype, copy=False)
    msg = f"unknown tropo mode {mode!r}"
    raise ValueError(msg)
