# ADR-0020: Keep the package layout; define the lean core by dependency tiers

- **Status:** accepted
- **Date:** 2026-10-06
- **Source:** plan task T17 (Venti package layout: lean core + extras)

## Decision

Venti keeps its existing packages (`gnss`, `spatial`, `filtering`, `solver`,
`unwrap`, `io`, `models`, `workflow`, `surface.py`). The "lean core" is
realised as dependency tiers in `pyproject.toml` — core, `[calibration]`,
`[decomposition]`, `[models]`, `[research]`, `[all]` — mirrored by pixi
features and the `ops` / `core-test` environments, not by moving modules into
a new `core/ calibration/ decomposition/` tree. Two structural changes are
made: the staging CLIs move from `scripts/staging/` into `venti.staging` (the
`sys.path` import hack goes away; thin wrappers stay), and
`models.plate_motion` / `models.load_itrf` are deleted in favour of
`geepers.euler`.

## Context

The plan's T17.1 proposed a reshuffled tree. Three facts argued against it:

1. cal-disp (the frozen operational SAS) imports `venti.surface`, `venti.gnss`,
   `venti.io`, `venti.workflow.config/utils` and `venti.log_setup` by name;
   a reshuffle would force a coordinated rename on the next pin bump for no
   functional gain.
2. The property the operational image needs — no dask/zarr/matplotlib/jupyter
   and no staging dependencies — is a property of what the modules import,
   not of where they live. The core modules already satisfy it once
   `opera-utils` is declared without its `disp` extra and the staging code
   is separated.
3. A 12 k-line move is unreviewable as a single-concern PR (CONTRIBUTING).

## Consequences

- `pip install venti` is the calibration engine; `venti[research]` adds
  staging, stack drivers' data access and notebooks (the former default).
- `tests/test_packaging.py` and the `core-only` CI job guard the tiers.
- The plan's T17.1/T17.2 wording is superseded by this record; later tasks
  (T18 `SensorSpec`, T28, T33) add modules beside the existing ones.
- `venti.gnss.unr` remains until T28; its duplication with geepers is known.
