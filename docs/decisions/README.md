# Decision records

Architecture decision records backfilled from the PRD decision log
(`docs/specs.md`, 2026-10-05). New decisions: copy `0000-template.md`.

| ADR | Title | Status |
|---|---|---|
| [ADR-0001](0001-venti-is-the-science-library-disp.md) | Venti is the science library, DISP-CAL the operational product | accepted |
| [ADR-0002](0002-disp-cal-delivers-a-per-pair.md) | DISP-CAL delivers a per-pair calibration surface; VLM is a separate product | accepted |
| [ADR-0003](0003-reference-frame-igs20-and-per-frame.md) | Reference frame IGS20 and per-frame plate lookup | accepted |
| [ADR-0004](0004-constant-velocity-unr-grid-is-operational.md) | Constant-velocity UNR grid is operational; variable grid is experimental | accepted |
| [ADR-0005](0005-the-interpolated-unr-grid-is-the.md) | The interpolated UNR grid is the calibration input; trust its sigma | accepted |
| [ADR-0006](0006-calibration-equals-the-exact-sum-of.md) | calibration equals the exact sum of its components | accepted |
| [ADR-0007](0007-unwrap-correction-is-gated-on-trade.md) | Unwrap correction is gated on trade study TS-U1 | accepted |
| [ADR-0008](0008-cli-and-runconfig-are-frozen-new.md) | CLI and runconfig are frozen; new layers come through output packaging | accepted |
| [ADR-0009](0009-defo-and-event-areas-are-curated.md) | Defo and event areas are curated GeoJSON; nodes and pixels inside are excluded | accepted |
| [ADR-0010](0010-geepers-owns-unr-access-euler-and.md) | geepers owns UNR access, Euler and GPS Imaging; tropo is shared; validation is separate | accepted |
| [ADR-0011](0011-extras-based-packaging-with-all.md) | Extras-based packaging with [all] | accepted |
| [ADR-0012](0012-sensor-agnostic-venti-s1-first-nisar.md) | Sensor-agnostic Venti; S1 first, NISAR as cal-disp v2 | accepted |
| [ADR-0013](0013-gamma-frozen-v0-5-calval-adopts.md) | Gamma frozen; v0.5 CalVal adopts the trade-study algorithm behind flags | accepted |
| [ADR-0014](0014-eight-frame-calval-benchmark.md) | Eight-frame CalVal benchmark | accepted |
| [ADR-0015](0015-vlm-design-per-pair-wls-with.md) | VLM design: per-pair, WLS with N from GNSS, provenance-checked E/N | accepted |
| [ADR-0016](0016-operations-one-run-per-granule-72.md) | Operations: one run per granule, 72 h latency, small EC2, frozen grid snapshot | accepted |
| [ADR-0017](0017-uncertainty-per-frame-k-sigma-disp.md) | Uncertainty: per-frame k, sigma_DISP from TS-S1, sigma realism reported not gated | accepted |
| [ADR-0018](0018-engineering-standards-claude-code-in-the.md) | Engineering standards, Claude Code in the loop, development on forks | accepted |
| [ADR-0019](0019-the-prd-lives-in-docs-specs.md) | The PRD lives in docs/specs.md | accepted |
| [ADR-0020](0020-keep-the-package-layout-tiers-by.md) | Keep the package layout; define the lean core by dependency tiers | accepted |
