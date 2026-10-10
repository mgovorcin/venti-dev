# Changelog

All notable changes to Venti. Format: [Keep a Changelog](https://keepachangelog.com/);
versions follow [PEP 440](https://peps.python.org/pep-0440/) via setuptools_scm.
Value-changing entries say which product layers move and by how much
(`workflow-regression` policy in `CONTRIBUTING.md`).

## [Unreleased]

### Changed

- **Dependency tiers (plan T17, ADR-0020).** `pip install venti` is now the
  calibration core (what cal-disp installs): `opera-utils` without its `disp`
  extra, so zarr/dask stay out of the operational image; staging, stack
  drivers' data access and notebooks are `venti[research]`;
  `venti[calibration]`, `[decomposition]` and `[models]` exist as tiers (empty
  today); `venti[all]` is everything. The former `notebooks` and `analysis`
  extras are folded into `research`. pixi: `default`/`test`/`docs`/`dev` keep
  the historical full set, `ops` is the core only, `core-test` core + tests.
- **Staging CLIs are a package.** `scripts/staging/{dem,disp,los,tropo,
  stage_frame,stage_window}_cli.py` and `utils.py` moved to `venti.staging`;
  the scripts remain as thin wrappers, so `python scripts/staging/dem_cli.py`
  still works and `venti run`'s staging no longer needs a source checkout
  (the `sys.path` import hack is gone). `scripts/staging/unr_cli.py` on
  upstream `main` is a stale copy of `los_cli.py` and was left untouched.
- **Plate motion comes from geepers.** `venti.models.plate_motion` and
  `venti.models.load_itrf` (and the ITRF JSON/CSV tables) are removed;
  `venti.models` re-exports `geepers.euler.plate_pole`, `plate_velocity_enu`,
  `load_plate_motion_model`, `PLATE_CODES`, `EulerPole`.
  `scripts/get_model_rates.py` uses them; its PMM raster attribute now says
  `mm/year`, which is what the values always were.
- `geepers[grid]` is a core dependency, pinned to a commit on the
  `mgovorcin/geepers` fork until it is tagged.
- pixi: `python >=3.11,<3.14` (several PyPI dependencies have no 3.14 wheels
  yet, and an unbounded floor let the solver pick 3.14); `geopandas-base`
  instead of `geopandas` on conda (the metapackage drags in matplotlib); and
  `opera-utils` comes from PyPI, because the conda-forge package has no extras
  and depends on asf_search/dask/zarr, which the `ops` environment must not
  carry.

- **Unwrap-error correction quantises region offsets in cycles of λ/2.**
  `estimate_calibration_surface(wavelength_m=...)` now takes the full radar
  wavelength (as written in the DISP product's `radar_wavelength`) and passes
  `wavelength_m / 2` to the corrector. Previously the full wavelength was used
  as the cycle, so corrections were twice the true cycle; this is why the
  gamma 0.3 delivery ships with `unwrap_error_correction: false`.
  **cal-disp:** `_unwrap_cycle_length_m()` halved the wavelength before
  calling Venti to work around this; when the Venti pin is bumped past this
  change, pass the full wavelength instead (plan T26/T37). No product values
  change while the flag is off.
- `UnwrapCorrector(cycle_length=...)` / `correct_region_offset(cycle_length=...)`
  replace the `wavelength` argument, whose name misdescribed the quantity.
  `wavelength=` still works for one release with a `DeprecationWarning` and
  keeps its numeric meaning (it is not halved).
- pixi manifest solves with pixi >= 0.48 (`dem-stitcher` and `geepers` were
  declared twice); `sqlite >= 3.45` floor so `import sqlite3` works in the
  solved environments; environments `dev` and `ops`; tasks `lint`, `golden`,
  `e2e`, `docs-serve`.

### Added

- `venti.calibration` package — the v0.5 engine next to the gamma `venti.surface`:
  - `gaps.fill_gaps(data, mask, method)` (plan T29, R-S3): nearest /
    nearest+smooth (default) / biharmonic / GDAL inverse-distance fill so the
    field is finite everywhere; `base_weights` gives filled pixels weight 0.02.
  - `loclin.loclin_surface(field, weights, pixel_m, cutoff_wavelength_meters)`
    (plan T30, R-S2): weighted local-linear Gaussian-kernel regression solved
    for every pixel from filtered moments; the kernel sigma follows from the
    half-response wavelength (`kernel_sigma_px`, response exactly 1/2 at the
    cutoff); low-coverage regions (sea, outside the swath) blend into a
    smoothed fallback so the surface is continuous and never 0.
  - `weights.fit_weights(field, base, sigma_px, coherence, coherence_power,
    robust)` (plan T31, R-S4): ``base x coherence^p x robust``; the robust
    factor down-weights pixels that disagree with their local low-pass level
    (Huber, or a 0/1 gate at 4 scaled MADs); a 10 cm blunder moves the 50 km
    surface by < 0.1 mm. With the gamma defaults it is the identity.
  - `remove_restore` (plan T32, A5/R-S5): `AreaDB`/`EventDB` from versioned
    GeoJSON (`defo_area_db_json`, `event_db_json`), events applied only to
    pairs that span `t0` (or overlap `[t0, t1]`), `remove_restore_mask` on the
    DISP grid, `sigma_inflation_inside` (sigma grows with distance into an
    area), `RemoveRestore.for_pair` with the versions and active event ids for
    the provenance. A 20 cm bowl inside an area leaves the surface flat there.
  - `two_pass.calibrate_pair(disp, gnss_los, mask, ref_point, options, pixel_m,
    cycle_m, ...)` (plan T33, R-S1/D6): the orchestrator. `surface.method =
    windowed_plane` reproduces the gamma estimator exactly (and keeps the
    unwrap shift as a component gamma lost); `loclin` runs gap fill → weights →
    local-linear kernel with optional two passes (robust frame-wide tie,
    unwrap hook on DISP − CAL₁, shifts applied, final surface). The result's
    `calibration` is exactly `cal_gnss_surface + cal_reference_offset +
    cal_tropo + cal_set + cal_unwrap_shift` (`assert_closed`).
  - `tropo` (plan T34, R-T1): the cal-disp tropo numerics
    (`interpolate_in_time`, `interpolate_to_dem_surface`,
    `compute_los_correction`, `pair_correction`) ported verbatim and kept
    bit-identical by a parity test that imports cal-disp when available;
    `stratified_tropo` / `stratified_delay` (per-epoch ``a + b x + c y + d h
    + e h**2`` fit, the height-dependent part only), `dem_relief` (p95 - p5),
    `choose_tropo_mode` (``legacy`` -> `apply_tropo_correction`; ``auto``:
    off below 1.5 km relief, stratified at or above; the frame table pins
    F16940 stratified and F08882/F08886/F08622 off) and `apply_tropo`, which
    returns the ``cal_tropo`` component for `calibrate_pair`. Schema v2
    defaults (`tropo.mode = legacy`) change nothing. `TropoOptions` now
    rejects `relief_off_meters > relief_stratified_meters`.
  - `uncertainty` (plan T35, R-E1/R-E2): `fit_sigma` (kernel-weighted
    residual variance over the kernel's effective pixel count
    ``4 pi sigma_px**2`` for unit weights), `sigma_cal` (``sqrt((k
    sigma_grid)**2 + sigma_fit**2 + sigma_tropo**2 + sigma_ref**2)`` times the
    inflation inside interpolated defo/event areas) and `resolve_k`
    (`uncertainty.k_grid`; `'frame_table'` must be applied first).
    `calibrate_pair` (loclin path) now returns `sigma_cal`, `sigma_fit` and
    `n_eff` in the units of `disp`, and takes `sigma_tropo` / `sigma_ref`
    terms (0 by default); the gamma path leaves them None. On synthetic data
    with white noise and a k-inflated grid error, (surface - truth) /
    sigma_CAL has std 0.8-1.25 (test). sigma_CAL covers the calibration
    only: users add DISP noise below the cutoff (R-E2).
- `venti.unwrap.regions` and `venti.unwrap.cycles` (plan T36, R-U1, D7): the
  v0.5 unwrap-error estimator next to the gamma `UnwrapCorrector`.
  `segment_regions(water_mask, valid_mask)` is the trade-study water-mask
  recipe (erode 3x, regrow, drop < 20 px, label connected land);
  `estimate_cycles(residual, labels, gnss_los, coherent, cycle_m, pixel_m,
  options)` visits regions largest first, measures each against anchored
  coherent land within 12 km on `DISP - CAL_1`, shifts whole cycles only
  when the jump is within 0.15 of an integer and the shift moves the region
  toward the GNSS field (veto), and lets only measured-0 or shifted regions
  anchor others; optional inversion-residual gate; free offsets refused.
  `apply_shifts` / `shift_field` give the `cal_unwrap_shift` component and
  `make_unwrap_hook` plugs it into `calibrate_pair`. `Decisions` records
  every region (CSV, summary). `UnwrapOptions` gains `anchor_distance_meters`,
  `cycle_tolerance`, `min_coherent_area_km2`, `min_edge_area_km2`,
  `residual_gate_cycles`. The 14-case F08882 bench (`tests/test_unwrap_bench.py`,
  needs `VENTI_UNWRAP_BENCH_DIR`) scores 14/14 with and without the gate.
  Still gated: `unwrap_error_correction` stays False until TS-U1.
- `venti.gnss.snapshot` + `scripts/snapshot_unr_grid.py` (plan T48.1, R-G4):
  frozen UNR grid snapshots `unr_grid_<version>_<date>/` with
  `snapshot.json`, `MANIFEST.sha256`, the lookup byte for byte and
  `nodes/<id>_<frame>_<grid_type>.tenv8` (the layout `GnssGridConfig`
  reads); retrieval through geepers with per-node retry and a failure list;
  `verify_snapshot`, `load_snapshot`, `grid_config_from_snapshot` (carries
  the snapshot id into the GNSS provenance). `docs/operations.md` documents
  the layout, storage (bucket TBD) and the roll-forward procedure (T48.4).
- `venti.decomposition` (plan T53, PRD section 3.2, R-E3): `decompose_wls(los,
  sigma, enu, north)` solves [E, U] per pixel from two or more looks with
  north fixed from GNSS, weights ``1 / (sigma_k^2 + n_k^2 sigma_N^2)``,
  closed-form 2x2 normal equations with sigma propagation and the
  condition number; `project_vertical(los, sigma, enu, east, north)` for one
  geometry with sigma propagation; `decompose(...)` dispatches per pixel
  (WLS where >= 2 looks and ``cond <= cond_max``, else an inverse-variance
  mean of the per-look projections) and returns `DecompositionResult` with
  the mode flag (none / wls / projection), `cond` and `n_looks`. The
  single-geometry validity assumption (long-wavelength horizontal motion
  only; unmodelled east leaks as ``e/u``) is documented in the module and
  tested. `workflow.decomposition.decompose_to_enu` now calls it (north
  fixed at 0, warning on ill-conditioned geometry) instead of returning
  zeros. The `[decomposition]` extra stays empty (numpy only).
- `venti.temporal.resample_timeseries(ts, dates, target_dates, window_days,
  method)` (plan T54, D15): per pixel, a robust linear fit to the epochs
  within the window evaluated at each target epoch (``huber_linear``
  default, ``theil_sen``, ``linear``), with the standard error of the line
  at the target from the window residuals; one-sided windows at the series
  ends, NaN beyond the span unless `max_extrapolation_days`; chunked over
  space. `resample_dataarray` wraps it for xarray (dask-parallel when the
  input is chunked). Tests cover denoising and sigma calibration, seasonal
  bias vs window length, the series ends, missing epochs, a blundered epoch
  and the xarray wrapper. The frame benchmark (T54.3) waits for the T56
  frame choice.
- `scripts/benchmark_golden_pair.py` + `docs/benchmarks.md` (plan T30.4): the
  F08882 golden pair calibrated with `calibrate_pair` for the gamma and the
  v0.5 option sets and scored at 123 real UNR stations
  (`tests/data/f08882_golden_pair_stations.csv`, from the trade study). Gamma
  600 km reproduces the trade-study 21.4 mm; loclin 50 km two-pass 13.6 mm
  (pair) / 12.3 mm (MIDAS); coherence/robust weights and the defo-area
  exclusion are within 1 mm on this metric; +2 mm bias vs gamma recorded for
  the e2e gate.
- `venti.gnss.sampling` (plan T28): `sample_gnss_enu(cfg, grid, exclude)` —
  the single GNSS sampling path for DISP-CAL and VLM. Nodes inside the frame
  bounds plus a buffer (R-G2), rates from the constant-grid tenv8 files
  (`variable` only in reprocessing mode, R-G1), optional exclusion of nodes
  inside defo/event areas with GPS-Imaging re-interpolation (R-G5), E/N/U and
  sigma fields on the target grid, and a `Provenance` record (grid version/
  type/frame, snapshot id, lookup hash, node counts, field hash, digest) for
  the product metadata (R-G4). `project_field_to_los` projects the gridded
  field per pixel with the LOS rasters, so nodes outside the swath need no
  LOS look (the trade-study LOS extrapolation is unnecessary).
- `venti.frames` (plan T27): versioned frame-parameter table
  (`src/venti/data/frame_parameters.json`) with the shape of cal-disp's
  `algorithm_parameters_overrides_json`; precedence default < frame <
  explicit; `materialized()` for the frozen cal-disp field;
  `calibration_options.frame` (plate validated against geepers' tables, name,
  benchmark category). Initial entries for F08882, F08886, F16940, F08622.
- `algorithm_parameters` schema v2 (plan T19): nested option groups
  `surface`, `weights`, `tropo`, `gnss`, `uncertainty`, `unwrap` under
  `calibration_options`, the three downsampling keys cal-disp already used,
  `schema_version` (a version-1 gamma file loads unchanged and is upgraded in
  memory), `extra = "forbid"` so typos fail at load time. Every new default
  reproduces the gamma behaviour. **`unwrap_error_correction` now defaults to
  off** (PRD R-U1; the gamma file sets it explicitly). Reference:
  `docs/algorithm_parameters.md`.
- `venti.sensor` (plan T18, PRD R-X1): `SensorSpec` with the radar wavelength
  and λ/2 cycle, product and static-layer filename grammars, mask/quality/
  correction layer names and readers (`read_wavelength` checks the product's
  `/identification/radar_wavelength` against the nominal value); `S1` is
  implemented, `NISAR` is registered and refused until cal-disp v2 (T58).
  `get_sensor(name)`, `sensor_for_file(path)`.
- SPDX license headers (BSD-3-Clause) on every Python file, enforced by the
  `spdx-header` pre-commit hook (`scripts/spdx_check.py`).
- Shared engineering kit: pre-commit (check-toml, nbstripout keeping outputs),
  GitHub issue/PR templates, `CONTRIBUTING.md`, `CLAUDE.md` working sections.
- mkdocs site (`mkdocs.yml`, `docs/`): PRD, implementation plan, task tracker,
  ADR-0001…0019, API reference; `docs.yaml` workflow deploys gh-pages.
- Tests: packaging regression (`tests/test_packaging.py`), cycle-length
  regression and surface integration (`tests/test_unwrap_corrections.py`,
  `tests/test_surface.py`), ITRF Euler-pole conversion
  (`tests/test_models_itrf.py`), cal-disp `algorithm_parameters.yaml`
  compatibility (`tests/test_config_compat.py`).

### Fixed

- `venti.models.load_itrf.convert_to_euler_poles` imported from a module path
  that did not exist (`plate_motion.euler_pole`) and raised `ImportError`.
- `calibrate_pair` passed `gnss_los_std` to the gamma windowed fit
  unconditionally; cal-disp passes it only with
  `weight_fit_by_gnss_uncertainty`. On the F08882 golden pair the gamma
  station RMSE moved from 21.4 mm (the trade-study value) to 24.3 mm with the
  sigma supplied; now gated like cal-disp (found by the T30.4 benchmark).
- The two-pass pass-1 tie is computed on a grid of at most 256 pixels across
  (`TIE_MAX_PX`) and interpolated back: its 4-frame-wide Gaussian moments on
  the 1288 x 1577 fit grid took about 4 of the 5 minutes of a full-frame run;
  the coarse and fine ties agree to well below a millimetre (tested).
- `sigma_cal` / `sigma_fit` / `n_eff` were factor-trimmed (7728 x 9462 for
  the 7733 x 9464 F08882 frame) and `calibrate_pair` failed at full
  resolution; each pixel now takes its block's node.

## [0.0.0] – upstream `2a7e61f` (2025-11)

The state of `opera-adt/Venti` `main` that cal-disp gamma 0.3 pins.
