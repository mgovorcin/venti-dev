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
- [ ] T13 geepers lean core — depends on T12
- [ ] T14 `geepers[grid]`: GPS Imaging + Euler, with exclusion areas — depends on T13
- [ ] T15 `geepers[analysis]` and `[all]` — depends on T13
- [x] T16 Venti bug fixes and packaging repair — done 2026-10-05 on `feature/venti-bugfixes` (T16.7 research-branch copy still open, see plan)
- [ ] T17 Venti package layout: lean core + extras — depends on T14, T16
- [ ] T18 `SensorSpec` abstraction — depends on T17
- [ ] T19 Venti algorithm-parameters schema — depends on T17
- [x] T20 Venti documentation site — done 2026-10-05 on `feature/docs-site` (gh-pages deploy verified only after first push to main)
- [ ] T21 Validation package: port the e2e core — depends on T11, T15
- [ ] T22 Station classes and per-class metrics — depends on T21
- [ ] T23 Gate logic and comparison modes — depends on T21
- [ ] T24 Reports and traceability — depends on T22, T23
- [ ] T25 Validation CLI, caching, batch execution — depends on T24
- [ ] T26 cal-disp foundations: pin Venti/geepers, drop duplicates, dependency budget — depends on T09, T13, T17
- [ ] T27 Frame-parameter table — depends on T19

## Phase 2: Science port into Venti

- [ ] T28 `sample_gnss_enu` — depends on T14, T18, T19
- [ ] T29 Gap filling and continuous surface support — depends on T19
- [ ] T30 Local-linear surface with physical cutoff — depends on T29
- [ ] T31 Robust coherence weights — depends on T30
- [ ] T32 Remove-restore: defo/event areas — depends on T29
- [ ] T33 Two-pass orchestration and component bookkeeping — depends on T30, T31, T32
- [ ] T34 Tropo modes — depends on T19, T27
- [ ] T35 σ_CAL model — depends on T27, T33
- [ ] T36 Unwrap-error module (gated) — depends on T18, T33
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
- [ ] T48 Frozen UNR grid snapshot — depends on T28
- [ ] T49 v0.5 Docker image, golden regeneration, changelog — depends on T46, T47, T48
- [ ] T50 VnV report and sign-off — depends on T24, T46
- [ ] T51 Release cal-disp v0.5 — depends on T49, T50

## Phase 5: VLM v0.1

- [ ] T52 VLM repo skeleton — depends on T07, T17
- [ ] T53 Venti `[decomposition]`: WLS and projection — depends on T18, T35
- [ ] T54 Temporal resampling of asc/desc — depends on T17
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
