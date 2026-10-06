# ADR-0010: geepers owns UNR access, Euler and GPS Imaging; tropo is shared; validation is separate

- **Status:** accepted
- **Date:** 2026-10-05
- **Source:** PRD decision log entry D10 (`docs/specs.md`)

## Decision

geepers is the single source for UNR access, plate motion (Euler) and GPS Imaging; tropo stays in cal-disp and is copied into Venti; separate validation package.

## Context

See `docs/specs.md`; the requirement IDs and trade studies that depend on this
decision are listed there. This record was backfilled from the PRD interview;
the discussion that led to it is summarised in the PRD sections it affects.

## Consequences

Changing this decision requires a new ADR that supersedes this one and an
update to `docs/specs.md`.
