# ADR-0005: The interpolated UNR grid is the calibration input; trust its sigma

- **Status:** accepted
- **Date:** 2026-10-05
- **Source:** PRD decision log entry D5 (`docs/specs.md`)

## Decision

The UNR interpolated grid is the calibration input; trust its σ; TS-G1 checks it against sites; no separate support layer.

## Context

See `docs/specs.md`; the requirement IDs and trade studies that depend on this
decision are listed there. This record was backfilled from the PRD interview;
the discussion that led to it is summarised in the PRD sections it affects.

## Consequences

Changing this decision requires a new ADR that supersedes this one and an
update to `docs/specs.md`.
