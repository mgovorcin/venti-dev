# Architecture

How Venti is laid out, which parts the operational products depend on, and
how the dependency tiers map onto that. Decisions that constrain this page are
[ADR-0010](decisions/0010-geepers-owns-unr-access-euler-and.md) (geepers is
the single source for GNSS access, Euler poles and GPS Imaging),
[ADR-0011](decisions/0011-extras-based-packaging-with-all.md) (extras-based
packaging) and [ADR-0020](decisions/0020-keep-the-package-layout-tiers-by.md)
(keep the package layout; tiers by dependency).

## Packages

```
venti/
  surface.py        estimate_calibration_surface(): the entry point cal-disp calls
  gnss/             reference grid -> station frames -> LOS (reference.py, los.py); unr.py (UNR
                    tenv8 readers; download helpers to be replaced by geepers in plan T28)
  spatial/          gap filling (GDAL), interpolation, resampling, SpatialProcessor
  filtering/        moving-window plane fit, Gaussian and low-pass filters
  solver/           design matrices, lscov, plane fitting
  temporal/         moving-window robust line fits that resample asc/desc time
                    series onto common epochs (T54)
  unwrap/           gamma corrector (unwrap_corrections) + v0.5 water-mask regions and
                    whole-cycle estimator (regions, cycles; gated, T36); cycle = lambda/2
  io/               rasters and NetCDF; CalProduct / VlmProduct writers
  models/           GIA rate rasters (ICE-6G_D, Caron 2018); plate motion re-exported from geepers
  workflow/         pydantic config (runconfig, algorithm parameters), stack-level calibration
                    and decomposition drivers, stage_frame_data
  staging/          frame-input staging CLIs (DISP-S1 download, DEM, LOS, tropo); thin wrappers
                    remain in scripts/staging/
  __main__.py       `venti config|run|run-single` (tyro)
```

### What the operational products import

cal-disp (gamma 0.3) imports exactly:

| import | role |
|---|---|
| `venti.surface.estimate_calibration_surface`, `SENTINEL1_WAVELENGTH_M` | the calibration surface |
| `venti.gnss.GNSSReference`, `compute_gnss_los`, `compute_gnss_los_std` | GNSS grid → LOS |
| `venti.io.read_netcdf_correction` | SET from the DISP product |
| `venti.workflow.config.CalibrationOptions` | algorithm options |
| `venti.workflow.utils.scratch_temp_dir`, `venti.log_setup.configure_logging` | utilities |

These names are stable; the plan's later work (T28 `sample_gnss_enu`, T33
two-pass `calibrate_pair`, T18 `SensorSpec`) adds entry points next to them
and cal-disp switches in T37.

## Dependency tiers

| Tier | Install | Contents | Added dependencies |
|---|---|---|---|
| **core** | `pip install venti` | everything above except `staging` and the stack drivers' data access | numpy, scipy, pandas, scikit-image, joblib, tqdm, xarray, rioxarray, netcdf4, h5netcdf, rasterio, gdal, geopandas, shapely, pyproj, requests, pydantic, pyyaml, tyro, opera-utils (no extras), rich, `geepers[grid]` |
| `[calibration]` | `venti[calibration]` | placeholder today (the core already is the calibration engine); numba etc. land here with the two-pass surface | — |
| `[decomposition]` | `venti[decomposition]` | `venti.decomposition`: per-pixel WLS for [E, U] with N from GNSS, single-geometry projection, mode flag (T53); numpy only, so the extra adds nothing today | — |
| `[models]` | `venti[models]` | GIA / PMM context rasters (`scripts/get_model_rates.py`) | — |
| `[research]` | `venti[research]` | `venti.staging`, stack drivers that need `opera_utils.disp`, notebooks | opera-utils[disp,asf] (zarr, dask, s3fs), asf_search, dem-stitcher, pyarrow, numba, matplotlib, ipykernel, ipywidgets, folium |
| `[all]` | `venti[all]` | all of the above | |

pixi mirrors the tiers: `ops` = core only (what the DISP-CAL image installs),
`core-test` = core + test tooling, `default`/`test`/`docs`/`dev` = core +
`research` (+ test/docs), i.e. the historical behaviour.

**Rule:** nothing in the core imports `opera_utils.disp` (zarr/dask), asf_search,
dem_stitcher, matplotlib or numba at module import time. `workflow.calibration`
imports `opera_utils.disp.rebase_reference` lazily inside the stack driver;
running a stack needs `[research]`. `tests/test_packaging.py` and the
`core-only` CI job check the contract.

## Why the layout was not reshuffled

The implementation plan (T17.1) sketched a `core/ calibration/ decomposition/
models/ workflow/` tree. The existing packages already separate those concerns,
cal-disp imports them by their current names, and the goal of T17 — a lean
core for the operational image — is a property of the dependency graph, not of
directory names. Moving 12 k lines would have made the PR unreviewable and
broken the pin without buying isolation. See ADR-0020.

## Upstream geepers as the single source

`venti.models.plate_motion` and `venti.models.load_itrf` were deleted in favour
of `geepers.euler` (`plate_pole`, `plate_velocity_enu`, the ITRF2020/2014
tables). `venti.models` re-exports them. The UNR download/reader helpers in
`venti.gnss.unr` still exist because `GNSSReference` and cal-disp's runconfig
depend on their on-disk layout; they are replaced by the geepers-backed
`sample_gnss_enu` in T28.
