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

## [0.0.0] – upstream `2a7e61f` (2025-11)

The state of `opera-adt/Venti` `main` that cal-disp gamma 0.3 pins.
