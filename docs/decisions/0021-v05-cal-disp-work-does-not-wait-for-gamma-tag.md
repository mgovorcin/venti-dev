# ADR-0021: v0.5 cal-disp work does not wait for the gamma tag

- **Status:** accepted
- **Date:** 2026-10-06
- **Source:** owner request to move on to CalVal and unblock the gated chain
  (T01 → T03–T06 → T09 → T26 → T37)

## Context

The plan made T09 (standards in cal-disp) depend on T06 (gamma released), and
T26/T37 hang off T09. T06 in turn waits on T01: Talib's sign-off on the gamma
configuration. As of 2026-10-06 the gamma release is **not** finished:
upstream PR opera-adt/cal-disp#21 is open and unreviewed, there is no `v0.3`
tag, the golden has not been rebuilt in Docker against an agreed
configuration, and the delivery documents and URLs are open (T03–T05). Every
remaining gamma step waits on someone outside this work.

The dependency existed only to keep the gamma release diff clean, not because
v0.5 code needs the released gamma.

## Decision

- v0.5 cal-disp work (T09, T26, T37) proceeds now on fork feature branches
  **stacked on `gamma-release`**, never on `gamma-release` itself and never
  in upstream PR #21.
- T06 still gates the **upstream** merge of that work: nothing v0.5 goes to
  opera-adt/cal-disp until gamma is tagged; the stack is rebased onto the
  tag when it exists.
- The gamma golden (`test_golden/`) stays the regression reference for every
  v0.5 branch with gamma-equivalent defaults (Venti invariant). If T01 changes
  the gamma configuration, the golden is rebuilt once and the stack rebased.
- T01.1 (the decision memo) is done: cal-disp `docs/decisions/0001-gamma-config.md`
  (fork PR mgovorcin/cal-disp#2), recommending 600 km, tropo off, unwrap off.

## Consequences

- Plan dependencies: T09 and T26 depend on T07/T13/T17 for development and
  on T06 only for the upstream PR.
- A gamma configuration change after T01.2 means one golden rebuild plus a
  rebase of the v0.5 stack; the cost is accepted.
