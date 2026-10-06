# ADR-0007: Unwrap correction is gated on trade study TS-U1

- **Status:** accepted
- **Date:** 2026-10-05
- **Source:** PRD decision log entry D7 (`docs/specs.md`)

## Decision

Unwrap correction gated on TS-U1 (islands with GNSS first, inversion residuals inland later).

## Context

See `docs/specs.md`; the requirement IDs and trade studies that depend on this
decision are listed there. This record was backfilled from the PRD interview;
the discussion that led to it is summarised in the PRD sections it affects.

## Consequences

Changing this decision requires a new ADR that supersedes this one and an
update to `docs/specs.md`.
