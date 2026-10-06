# Development notes

Baseline recorded 2026-10-05 on `feature/venti-bugfixes` (plan T16.6). Update
the numbers when they change on purpose.

## Environments

| pixi env | Purpose | Installed here |
|---|---|---|
| `default` | Venti + runtime deps | yes |
| `dev` | default + test/lint/docs tooling (`test` + `docs` features) | yes |
| `test`, `docs`, `analysis` | upstream's narrower envs, kept for compatibility | yes / yes / no |
| `ops` | the default dependency set only — what an operational image installs | lock only |

```bash
export RATTLER_CACHE_DIR=/u/aurora-r0/govorcin/.cache/rattler UV_CACHE_DIR=/u/aurora-r0/govorcin/.cache/uv
export TMPDIR=/u/aurora-r0/govorcin/tmp
pixi install -e dev
pixi run -e dev test                       # pytest with coverage + doctests
pixi run -e dev lint                       # pre-commit run -a
JUPYTER_PREFER_ENV_PATH=1 pixi run -e docs docs   # mkdocs build --strict
```

## Baseline

| Check | Result (2026-10-05) |
|---|---|
| `pytest` (upstream `2a7e61f`, dev env) | 257 passed; coverage 74 % |
| `pytest` after T08/T16 | see CHANGELOG; all new tests pass, no upstream test changed its outcome |
| `pre-commit run -a` | clean (ruff, black, mypy, nbstripout, SPDX) |
| `mkdocs build --strict` | clean |

## Traps on this machine

- **pixi 0.48 rejects duplicate declarations** of a package in
  `[project.dependencies]` and `[tool.pixi.*pypi-dependencies]`;
  `tests/test_packaging.py` guards it.
- **sqlite clobber:** without a floor the solver picked `sqlite 3.32.3`, which
  overwrote `libsqlite3.so` and broke `import sqlite3`. Pinned `>= 3.45`.
- **PyPI timeouts** are frequent; `pixi install` may need two or three runs.
- **nbconvert templates:** stale templates under `~/.local/share/jupyter`
  shadow the env's and make mkdocs-jupyter fail with
  `'mermaid_js' is undefined`; set `JUPYTER_PREFER_ENV_PATH=1`.
- **Working directory:** run `pytest` / `pre-commit` from the repo root;
  `pixi run --manifest-path ... <cmd>` does not change the cwd.
- Full-frame runs need `NUMPY_MADVISE_HUGEPAGE=0`.

## Packaging findings (T16.1)

- Every third-party import in `src/venti` is declared in `pyproject.toml`
  (checked with an `ast` scan; aliases `osgeo`→`gdal`, `yaml`→`pyyaml`,
  `skimage`→`scikit-image`).
- `src/venti/workflow/stage_frame_data.py` imports `dem_cli`, `disp_cli`,
  `los_cli`, `tropo_cli` from `scripts/staging/` through a `sys.path` insert.
  Those modules are not part of the installed package, so the stage workflow
  only works from a source checkout. **Resolved in T17:** the CLIs are
  `venti.staging` (shims remain in `scripts/staging/`); their dependencies are
  the `research` extra. `scripts/staging/unr_cli.py` on upstream `main` is a
  stale copy of `los_cli.py` (same functions, one comment differs) and was
  left in place; it is not imported anywhere.
- `README.md` referenced `environment.yml`, which exists; the pixi route is
  now documented alongside it.

## Coordination with cal-disp

cal-disp passes `wavelength_m=_unwrap_cycle_length_m(λ)` (λ/2) to
`estimate_calibration_surface` to work around the former λ-vs-λ/2 bug. After
the fix in this branch Venti halves internally, so cal-disp must pass the
full wavelength when it bumps its Venti pin (plan T26/T37). With unwrap
correction off (gamma 0.3), product values are unaffected either way.
