# TODO tracker

Task-level status for [plan.md](plan.md). One line per task; subtask checkboxes
live in `plan.md`. A task is checked here only when every subtask in `plan.md`
is checked and the definition of done in `plan.md` holds. "Blocked" names the
external input that is missing.

Legend: `[x]` done · `[~]` in progress · `[ ]` not started · `[!]` blocked

## Phase 0: Gamma 0.3 release

- [!] T01 Confirm the gamma configuration with Talib — blocked on Talib's answers
- [x] T02 Secure a Docker build host and build the gamma image — done 2026-10-05: aurora has Docker; `cal-disp:0.3.0-rc` built, record on cal-disp `feature/docker-build-record`
- [ ] T03 Rebuild the golden inside Docker and validate — depends on T01, T02
- [ ] T04 Fix the delivery documents — depends on T01
- [ ] T05 Confirm product URLs and version string — depends on T01
- [ ] T06 Release gamma 0.3 — depends on T03, T04, T05

## Phase 1: Foundations

- [x] T07 Shared engineering-standards kit
- [x] T08 Apply standards to Venti (`venti-dev`) — done 2026-10-05 on `feature/standards` (pre-commit.ci switch-on is a manual owner step)
- [ ] T09 Apply standards to cal-disp — depends on T06, T07
- [x] T10 Apply standards to geepers fork — done 2026-10-05 on geepers `feature/standards` (3 commits: env, kit, type fixes)
- [ ] T11 Create the validation package repo — depends on T07
- [~] T12 geepers dependency audit — T12.1–T12.2 done (`docs/dependency_audit.md` on `feature/extras-split`); T12.3 owner sign-off pending
- [x] T13 geepers lean core — done 2026-10-06 on geepers `feature/extras-split` (provisional on the T12.3 sign-off)
- [x] T14 `geepers[grid]`: GPS Imaging + Euler, with exclusion areas — done 2026-10-06 on geepers `feature/extras-split` (reinterpolate_nodes, plate tables, plate_velocity_enu)
- [x] T15 `geepers[analysis]` and `[all]` — done 2026-10-06 on geepers `feature/extras-split` (CI matrix per tier, README install matrix, CHANGELOG)
- [x] T16 Venti bug fixes and packaging repair — done 2026-10-05 on `feature/venti-bugfixes` (T16.7 research-branch copy still open, see plan)
- [x] T17 Venti package layout: lean core + extras — done 2026-10-06 on `feature/lean-core` (tiers by dependency, ADR-0020; `venti.staging`; plate motion from geepers; `gnss/unr.py` deferred to T28)
- [x] T18 `SensorSpec` abstraction — done 2026-10-06 on `feature/sensor-spec` (`venti.sensor`: S1 implemented, NISAR registered/refused until T58)
- [x] T19 Venti algorithm-parameters schema — done 2026-10-06 on `feature/algorithm-schema` (schema v2: nested option groups, versioned loader, extra=forbid, unwrap default off)
- [x] T20 Venti documentation site — done 2026-10-05 on `feature/docs-site` (gh-pages deploy verified only after first push to main)
- [ ] T21 Validation package: port the e2e core — depends on T11, T15
- [ ] T22 Station classes and per-class metrics — depends on T21
- [ ] T23 Gate logic and comparison modes — depends on T21
- [ ] T24 Reports and traceability — depends on T22, T23
- [ ] T25 Validation CLI, caching, batch execution — depends on T24
- [ ] T26 cal-disp foundations: pin Venti/geepers, drop duplicates, dependency budget — depends on T09, T13, T17
- [x] T27 Frame-parameter table — done 2026-10-06 on `feature/frame-table` (`venti.frames`, `calibration_options.frame`, bundled table for 4 frames; T27.4 cal-disp consumption waits for T37)

## Phase 2: Science port into Venti

- [x] T28 `sample_gnss_enu` — done 2026-10-06 on `feature/gnss-sampling` (`venti.gnss.sampling`: buffer, exclusion + re-interpolation, E/N/U fields, provenance; per-pixel LOS projection replaces LOS extrapolation)
- [x] T29 Gap filling and continuous surface support — done 2026-10-06 on `feature/calibration-surface` (`venti.calibration.gaps`)
- [~] T30 Local-linear surface with physical cutoff — done except T30.4 (golden-pair benchmark after the cal-disp wiring, T37); method switch in `calibrate_pair`
- [~] T31 Robust coherence weights — core done 2026-10-06 (`venti.calibration.weights`); gamma-mask golden regression (T31.3) is a cal-disp run in T37
- [x] T32 Remove-restore: defo/event areas — done 2026-10-06 on `feature/calibration-weights-rr` (`venti.calibration.remove_restore`)
- [x] T33 Two-pass orchestration and component bookkeeping — done 2026-10-06 on `feature/calibration-two-pass` (`venti.calibration.two_pass.calibrate_pair`, `CalibrationResult`; gamma path reproduces `estimate_calibration_surface` to 1e-7)
- [x] T34 Tropo modes — done 2026-10-06 on `feature/tropo-modes` (`venti.calibration.tropo`: cal-disp numerics bit-identical, stratified fit, relief rule, `apply_tropo` → `cal_tropo`; cal-disp CI import deferred to T37)
- [x] T35 σ_CAL model — done 2026-10-06 on `feature/sigma-cal` (`venti.calibration.uncertainty`; `calibrate_pair` returns `sigma_cal`; z-score test 0.8–1.25)
- [x] T36 Unwrap-error module (gated) — done 2026-10-06 on `feature/unwrap-module` (`venti.unwrap.regions/cycles`; 14/14 bench; still off by default until TS-U1)
- [ ] T37 cal-disp wiring to the new Venti workflow — depends on T26, T33, T34, T35
- [ ] T38 e2e on the 4 existing frames — depends on T25, T37

## Phase 3: Trade studies and benchmark data

- [ ] T39 TS-G1: grid fidelity and per-frame k — depends on T21, T28
- [ ] T40 Benchmark data staging (frames 5–8) — depends on T11
- [ ] T41 Curate defo and event GeoJSON databases — depends on T32, T40
- [ ] T42 TS-U1 phase 1 — depends on T36, T38
- [ ] T43 TS-S1: DISP noise model — depends on T21, T38
- [ ] T44 TS-T1 and TS-B1 — depends on T31, T34, T38

## Phase 4: cal-disp v0.5 CalVal release

- [ ] T45 Populate the frame table — depends on T27, T39, T44
- [ ] T46 8-frame gate and edge cases — depends on T38, T40, T41, T45
- [ ] T47 Memory and runtime toward a small EC2 instance — depends on T37
- [~] T48 Frozen UNR grid snapshot — T48.1 + T48.4 done 2026-10-06 on `feature/unr-snapshot` (`venti.gnss.snapshot`, `scripts/snapshot_unr_grid.py`, `docs/operations.md`); T48.2 blocked on the S3 bucket decision, T48.3 waits for T37
- [ ] T49 v0.5 Docker image, golden regeneration, changelog — depends on T46, T47, T48
- [ ] T50 VnV report and sign-off — depends on T24, T46
- [ ] T51 Release cal-disp v0.5 — depends on T49, T50

## Phase 5: VLM v0.1

- [ ] T52 VLM repo skeleton — depends on T07, T17
- [x] T53 Venti `[decomposition]`: WLS and projection — done 2026-10-06 on `feature/decomposition` (`venti.decomposition`: `decompose_wls`, `project_vertical`, `decompose` with mode flag)
- [~] T54 Temporal resampling of asc/desc — T54.1–T54.2 done 2026-10-06 on `feature/temporal-resampling` (`venti.temporal`); T54.3 benchmark waits for the T56 frame choice
- [ ] T55 GNSS E/N for VLM with provenance check — depends on T28, T37
- [ ] T56 VLM workflow and product writer — depends on T53, T54, T55
- [ ] T57 VLM validation and v0.1 release — depends on T24, T56

## Phase 6: DISP-NISAR

- [ ] T58 NISAR `SensorSpec` — depends on T18, beta data
- [ ] T59 TS-N1: ionosphere — depends on T38, T58
- [ ] T60 NISAR benchmark and cal-disp v2 — depends on T59

## Later

- [ ] T61 Output packaging: component and GNSS layers — depends on T14, T51
- [ ] T62 TS-U1 phase 2: inland unwrap errors — depends on T42

## Log

- 2026-10-05 Tracker created. Phase 0 is blocked on external inputs (Talib, a Docker host); starting phase 1 at T07.
- 2026-10-05 T07 done (`00_tools/standards`, 16 kit tests). T08 done: Venti `feature/prd-docs` → `feature/standards` (pixi fix, kit, SPDX) → `feature/docs-site` (T20) and `feature/venti-bugfixes` (T16: λ/2 cycle, ITRF import, cal-disp YAML compat; 271 tests, lint clean). T10/T12 in progress on geepers `feature/standards` (dependency audit written; `affine` 3 warning filter + mypy hook bump needed for a green dev env). T02 unblocked: Docker daemon available on aurora; gamma image build started (`cal-disp:0.3.0-rc`).
- 2026-10-05 (cont.) geepers `feature/standards` committed (env, kit, 21 type fixes; 348 tests, lint green); `feature/extras-split` holds the T12 audit. Venti branches rebased into one linear stack (prd-docs → standards → docs-site → venti-bugfixes). Gamma image `cal-disp:0.3.0-rc` built and validated in-image at 1e-6 (cal-disp `feature/docker-build-record`). Open decisions: T12.3 partition sign-off; validation/VLM repo names (T11/T52); T01 with Talib.
- 2026-10-06 All branches pushed; PRs open on the forks (venti-dev #1–#4, geepers #1–#2, cal-disp #1). T13 done: geepers core imports without geopandas/pandera, deps split core/[grid]/[analysis]/[plot]/[all], `core-only` CI job. Fork Actions are off until enabled in each fork's Actions tab. Next: T14 (incl. T14.3a ITRF2020 plate table), T15.
- 2026-10-06 T14 and T15 done on geepers `feature/extras-split` (PR #2): ITRF2020/2014 plate tables + `plate_velocity_enu`, `reinterpolate_nodes`, CI matrix per tier, README/CHANGELOG. 376 tests. Phase 1 remaining: T09 (blocked on T06), T11 (repo name), T17–T19, T21–T27. Next: T17 (Venti layout) now that T14 and T16 are done.
- 2026-10-06 T17 done on `feature/lean-core` (PR #5): tiers by dependency (ADR-0020), `venti.staging`, plate motion from geepers, `ops`/`core-test` envs; suite 291 passed in 19 s. Next: T18 `SensorSpec`, T19 algorithm-parameters schema, T27 frame table.
