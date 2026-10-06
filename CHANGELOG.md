# Changelog

All notable changes to Venti. Format: [Keep a Changelog](https://keepachangelog.com/);
versions follow [PEP 440](https://peps.python.org/pep-0440/) via setuptools_scm.
Value-changing entries say which product layers move and by how much
(`workflow-regression` policy in `CONTRIBUTING.md`).

## [Unreleased]

### Changed

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
