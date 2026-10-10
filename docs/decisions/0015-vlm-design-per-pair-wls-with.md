# ADR-0015: VLM design: per-pair, WLS with N from GNSS, provenance-checked E/N

- **Status:** accepted
- **Date:** 2026-10-05
- **Source:** PRD decision log entry D15 (`docs/specs.md`)

## Decision

VLM: per-pair; asc+desc WLS with N from GNSS, projection otherwise; moving-window temporal resampling; E/N from DISP-CAL or the grid, with a provenance check.

## Context

See `docs/specs.md`; the requirement IDs and trade studies that depend on this
decision are listed there. This record was backfilled from the PRD interview;
the discussion that led to it is summarised in the PRD sections it affects.

## Consequences

Changing this decision requires a new ADR that supersedes this one and an
update to `docs/specs.md`.
