# Python Coding & Editing Guidelines

> **Living document – PRs welcome!**
> Last updated: 2025‑07‑15

## Table of Contents

1. Philosophy
1. Docstrings & Comments
1. Type Hints
1. Documentation

---

## Philosophy

- **Readability, reproducibility, performance – in that order.**
- Prefer explicit over implicit; avoid hidden state and global flags.
- Measure before you optimize (`time.perf_counter`, `line_profiler`).
- Each module holds a **single responsibility**; keep public APIs minimal.

## Docstrings & Comments

- Style: NumPyDoc.
- Start with a one‑sentence summary in the imperative mood.
- Sections: Parameters, Returns, Raises, Examples, References.
- Use backticks for code or referring to variables (e.g. `xarray.DataArray`).
- Do not use emojis, or non-unicode characters in comments/print statements.
- Cite peer‑reviewed papers with DOI links when relevant.
- Write code that explains itself rather than needs comments.
- For the inline you do add, explain *why*, not what. For example, *don't* write:

```python
# open the file
f = open(filename)
```

- The comments should be things which are not obvious to a reader with typical background knowledge.

## Tools

- ruff is use for most code maintenance, black for formatting, mypy for type checking, pytest for testing
- You can run `pre-commit run -a` to run all pre-commit hooks and check for style violations

## Code Style

- Annotate all public functions (PEP 484).
- Prefer `Protocol` over `ABC`s when only an interface is needed.
- Validate external inputs via Pydantic models (if existing); otherwise, use `dataclasses`
- Parse, don't validate, with your dataclasses. Checks should be at the serialization boundaries, not scattered everywhere in the code.
- If you need to add an ignore, ignore a specific check like # type: ignore[specific]
- Don't write error handing code or smooth over exceptions/errors unless they are expected as part of control flow.
- In general, write code that will raise an exception early if something isn't expected.
- Enforce important expectations with asserts, but raise errors for user-facing problems.

## Documentation

- mkdocs + Jupyter. Hosted on ReadTheDocs.
- Auto API from type hints.
- Provide tutorial notebooks covering common workflows.
- Include examples in docstrings.
- Add high-level guides for key functionality.

---

# Working in this repo

Sections below follow `00_tools/standards/CLAUDE.md.template`. The PRD is
`docs/specs.md`; the task plan is `docs/plan.md`; status is `docs/todo.md`.

## What this repo is

Venti is the sensor-agnostic InSAR–GNSS fusion library: GNSS grid → LOS
sampling, the calibration surface, unwrap-cycle correction, tropo modes,
remove-restore, and (later) LOS decomposition / projection to vertical. The
operational products (`cal-disp`, VLM) call it; Venti itself owns no product
format, runconfig or delivery packaging.

## Architecture (current `main`)

`gnss/` reference grid and LOS projection · `spatial/` gap filling,
interpolation, resampling · `filtering/` moving-window plane · `surface.py`
the calibration surface entry point used by `cal-disp` · `unwrap/` region
cycle correction · `models/` ITRF PMM and GIA rates · `io/` rasters · `workflow/`
stack-level drivers · `__main__.py` the `venti` tyro CLI. The target layout
(core + `[calibration]`, `[decomposition]`, `[models]`, `[research]` extras) is
in `docs/plan.md` T17.

## Invariants (do not break without a decision record in `docs/decisions/`)

- Sensor parameters (wavelength, readers) come from a `SensorSpec`; never
  hard-code Sentinel-1 values in algorithm code. One LOS displacement cycle is
  **λ/2**, not λ.
- `calibration == Σ cal_* components` to 1e-6 (once components exist, T33).
- Science changes land behind an `algorithm_parameters` flag whose default
  reproduces the previous result; `cal-disp`'s golden must stay green until a
  deliberate regeneration.
- No heavy dependencies (dask, zarr, matplotlib, jupyter) in the core
  dependency set; they belong to extras.

## Commands

```bash
export RATTLER_CACHE_DIR=/u/aurora-r0/govorcin/.cache/rattler UV_CACHE_DIR=/u/aurora-r0/govorcin/.cache/uv
export TMPDIR=/u/aurora-r0/govorcin/tmp
pixi install -e dev
pixi run -e dev test      # pytest (doctests included)
pixi run -e dev lint      # pre-commit run -a (ruff, black, mypy, nbstripout, SPDX)
JUPYTER_PREFER_ENV_PATH=1 pixi run -e docs docs   # mkdocs build --strict; the env var stops stale
                                                 # ~/.local nbconvert templates from shadowing the env's
```

Full-frame runs on aurora need `NUMPY_MADVISE_HUGEPAGE=0`.

## Golden / value-change policy

Venti has no golden of its own; `cal-disp` does. A Venti change that can alter
`cal-disp` output is tested there (`pixi run validate --golden-dir test_golden`
in `00_tools/src/cal-disp`) before the Venti PR is opened. See the
`workflow-regression` skill.

## Tests

Every change ships unit tests in the same commit. Numeric functions also get a
synthetic-truth test (known plane / bowl / cycle jump recovered within
tolerance), and anything that shapes a run gets a regression test pinning its
output. See the `unit-tests` skill and `CONTRIBUTING.md`.

## Conventions

- Branches: never commit on `main`; `feature/<topic>` for one concern; PRs to
  `mgovorcin/venti-dev` `main`, reviewed by the owner; upstream PRs to
  `opera-adt/Venti` only at milestones.
- SPDX header (`BSD-3-Clause`) on every `.py` file; `scripts/spdx_check.py`
  enforces it via pre-commit.
