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
