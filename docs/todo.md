# TODO tracker

Task-level status for [plan.md](plan.md). One line per task; subtask checkboxes
live in `plan.md`. A task is checked here only when every subtask in `plan.md`
is checked and the definition of done in `plan.md` holds. "Blocked" names the
external input that is missing.

Legend: `[x]` done · `[~]` in progress · `[ ]` not started · `[!]` blocked

## Phase 0: Gamma 0.3 release

- [x] T01 Confirm the gamma configuration with Talib — done 2026-10-06: 600 km, tropo off, unwrap off (cal-disp `docs/decisions/0001-gamma-config.md`)
- [x] T02 Secure a Docker build host and build the gamma image — done 2026-10-05: aurora has Docker; `cal-disp:0.3.0-rc` built, record on cal-disp `feature/docker-build-record`
- [x] T03 Rebuild the golden inside Docker and validate — done 2026-10-06: `cal-disp:0.3.0-rc2`, 1e-6 in image and on aurora, manifest committed (`8bafac7`)
- [ ] T04 Fix the delivery documents — depends on T01
- [ ] T05 Confirm product URLs and version string — depends on T01
- [ ] T06 Release gamma 0.3 — depends on T03, T04, T05

## Phase 1: Foundations

- [x] T07 Shared engineering-standards kit
- [x] T08 Apply standards to Venti (`venti-dev`) — done 2026-10-05 on `feature/standards` (pre-commit.ci switch-on is a manual owner step)
- [ ] T09 Apply standards to cal-disp — depends on T07; T06 only for the upstream PR (ADR-0021)
- [x] T10 Apply standards to geepers fork — done 2026-10-05 on geepers `feature/standards` (3 commits: env, kit, type fixes)
- [x] T11 Create the validation package repo — done 2026-10-06: private `mgovorcin/disp2vlm_validation` (`v0.0.0`)
- [x] T12 geepers dependency audit — done 2026-10-06 (T12.3 approved by the owner on geepers PR #2)
- [x] T13 geepers lean core — done 2026-10-06 on geepers `feature/extras-split` (provisional on the T12.3 sign-off)
- [x] T14 `geepers[grid]`: GPS Imaging + Euler, with exclusion areas — done 2026-10-06 on geepers `feature/extras-split` (reinterpolate_nodes, plate tables, plate_velocity_enu)
- [x] T15 `geepers[analysis]` and `[all]` — done 2026-10-06 on geepers `feature/extras-split` (CI matrix per tier, README install matrix, CHANGELOG)
- [x] T16 Venti bug fixes and packaging repair — done 2026-10-05 on `feature/venti-bugfixes` (T16.7 research-branch copy still open, see plan)
- [x] T17 Venti package layout: lean core + extras — done 2026-10-06 on `feature/lean-core` (tiers by dependency, ADR-0020; `venti.staging`; plate motion from geepers; `gnss/unr.py` deferred to T28)
- [x] T18 `SensorSpec` abstraction — done 2026-10-06 on `feature/sensor-spec` (`venti.sensor`: S1 implemented, NISAR registered/refused until T58)
- [x] T19 Venti algorithm-parameters schema — done 2026-10-06 on `feature/algorithm-schema` (schema v2: nested option groups, versioned loader, extra=forbid, unwrap default off)
- [x] T20 Venti documentation site — done 2026-10-05 on `feature/docs-site` (gh-pages deploy verified only after first push to main)
- [x] T21 Validation package: port the e2e core — depends on T11, T15
- [x] T22 Station classes and per-class metrics — depends on T21
- [x] T23 Gate logic and comparison modes — depends on T21
- [x] T24 Reports and traceability — depends on T22, T23
- [x] T25 Validation CLI, caching, batch execution — depends on T24 (T25.3 batch deferred to T46)
- [~] T26 cal-disp foundations: pin Venti/geepers, drop duplicates, dependency budget — T26.1/T26.3 (wavelength) PR #3, T26.6 PR #4, T26.4–T26.5 PR #5 done; T26.2 swap waits for the next pin bump (geepers hardened in PR #4); T26.3 SensorSpec wrappers open
- [x] T27 Frame-parameter table — done 2026-10-06 on `feature/frame-table` (`venti.frames`, `calibration_options.frame`, bundled table for 4 frames; T27.4 cal-disp consumption waits for T37)

## Phase 2: Science port into Venti

- [x] T28 `sample_gnss_enu` — done 2026-10-06 on `feature/gnss-sampling` (`venti.gnss.sampling`: buffer, exclusion + re-interpolation, E/N/U fields, provenance; per-pixel LOS projection replaces LOS extrapolation)
- [x] T29 Gap filling and continuous surface support — done 2026-10-06 on `feature/calibration-surface` (`venti.calibration.gaps`)
- [x] T30 Local-linear surface with physical cutoff — done 2026-10-06; T30.4 benchmark on `feature/golden-benchmark` (`docs/benchmarks.md`: 21.4 → 13.6 mm station RMSE on the F08882 golden pair)
- [~] T31 Robust coherence weights — core done 2026-10-06 (`venti.calibration.weights`); gamma-mask golden regression (T31.3) is a cal-disp run in T37
- [x] T32 Remove-restore: defo/event areas — done 2026-10-06 on `feature/calibration-weights-rr` (`venti.calibration.remove_restore`)
- [x] T33 Two-pass orchestration and component bookkeeping — done 2026-10-06 on `feature/calibration-two-pass` (`venti.calibration.two_pass.calibrate_pair`, `CalibrationResult`; gamma path reproduces `estimate_calibration_surface` to 1e-7)
- [x] T34 Tropo modes — done 2026-10-06 on `feature/tropo-modes` (`venti.calibration.tropo`: cal-disp numerics bit-identical, stratified fit, relief rule, `apply_tropo` → `cal_tropo`; cal-disp CI import deferred to T37)
- [x] T35 σ_CAL model — done 2026-10-06 on `feature/sigma-cal` (`venti.calibration.uncertainty`; `calibrate_pair` returns `sigma_cal`; z-score test 0.8–1.25)
- [x] T36 Unwrap-error module (gated) — done 2026-10-06 on `feature/unwrap-module` (`venti.unwrap.regions/cycles`; 14/14 bench; still off by default until TS-U1)
- [~] T37 cal-disp wiring to the new Venti workflow — T37.1, T37.2, T37.4–T37.6 done 2026-10-09 (cal-disp PR #6; golden 1e-6, v0.5 closure exact, 55 s / 7.27 GB); T37.3 metadata needs a deliberate golden update
- [x] T38 e2e on the 4 existing frames — depends on T25, T37

## Phase 3: Trade studies and benchmark data

- [x] T39 TS-G1: grid fidelity and per-frame k — depends on T21, T28 (σ₀ term proposal open)
- [ ] T40 Benchmark data staging (frames 5–8) — depends on T11
- [ ] T41 Curate defo and event GeoJSON databases — depends on T32, T40
- [x] T42 TS-U1 phase 1 — depends on T36, T38 (no-go for v0.5; 3 bugs fixed on the way)
- [x] T43 TS-S1: DISP noise model — depends on T21, T38
- [~] T44 TS-T1 and TS-B1 — TS-B1 done 2026-10-09 (keep p = 8); TS-T1 waits for a 0.3–1.5 km relief frame (T40)

## Phase 4: cal-disp v0.5 CalVal release

- [ ] T45 Populate the frame table — depends on T27, T39, T44
- [ ] T46 8-frame gate and edge cases — depends on T38, T40, T41, T45
- [~] T47 Memory and runtime toward a small EC2 instance — v0.5 7.27 → 4.1 GB, budget 4.5 GB (cal-disp PR #8, venti-dev PR #22); EC2 measurement open
- [~] T48 Frozen UNR grid snapshot — T48.1, T48.2, T48.4 done 2026-10-06: `unr_grid_0.3_IGS20_20261006` on `s3://opera-adt/opera-ancillary/unr-grid/` (28,492 nodes); T48.3 waits for T37
- [ ] T49 v0.5 Docker image, golden regeneration, changelog — depends on T46, T47, T48
- [ ] T50 VnV report and sign-off — depends on T24, T46
- [ ] T51 Release cal-disp v0.5 — depends on T49, T50

## Phase 5: VLM v0.1

- [~] T52 VLM repo skeleton — T52.1 done 2026-10-06: private `mgovorcin/opera_vlm` (VLM-S1 + VLM-NI, scaffold `v0.0.0`, PR #1 invariants); T52.2 runconfig ADR and T52.3 product spec open
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
- 2026-10-06 Phase 2 science port finished on the Venti side: T34 tropo modes (PR #13), T35 σ_CAL (PR #14), T36 unwrap module, gated, 14/14 bench (PR #15); then the unblocked later-phase work: T48.1/T48.4 UNR snapshot + operations doc (PR #16), T53 decomposition (PR #17), T54.1–T54.2 temporal resampling (PR #18). Suite 456 passed, 7 skipped (parity/golden/bench tests need `PYTHONPATH`, `CAL_DISP_GOLDEN_DIR`, `VENTI_UNWRAP_BENCH_DIR`). Everything left open needs a decision or an upstream step: T12.3 sign-off; T11/T52 repo names; T01 → T03–T06 → T09 → T26 → T37 (cal-disp wiring, which also closes T27.4, T30.4, T31.3, T34 CI import, T48.3); T48.2 S3 bucket; T58 beta DISP-NISAR data; fork Actions still to be enabled.
- 2026-10-06 (later) CI on the fork fixed (full matrix installs every tier; pre-commit under 3.13; files committed before linting formatted at every branch tip). T30.4 benchmark done on `feature/golden-benchmark` (PR #19): gamma path reproduces the trade study exactly (21.4 mm); loclin 50 km 13.6 mm; three bugs found by the full-frame run and fixed mid-stack (gamma GNSS-sigma weights gate, sigma-map trimming, 290 s pass-1 tie). Open from it: TS-B1 on weights, +2 mm bias for T38/T39.
- 2026-10-06 Gamma release confirmed **not finished** (upstream PR #21 open and unreviewed, no v0.3 tag; T01 sign-off, T03 Docker golden, T04 docs, T05 URLs open). T01.1 memo written (cal-disp fork PR #2: recommends 600 km, tropo off, unwrap off). ADR-0021: v0.5 cal-disp work proceeds on fork branches stacked on `gamma-release`; T06 gates only the upstream merge. T11 done (`disp2vlm_validation`, private). T12.3 sign-off requested on geepers PR #2. T48.2 proposal: `s3://opera-adt/opera-ancillary/unr-grid/<snapshot_id>/`.
- 2026-10-06 (owner decisions) Talib accepted the gamma config (600 km, tropo off, unwrap off) → T01 done; golden rebuilt and validated at 1e-6 inside `cal-disp:0.3.0-rc2` and on aurora, manifest committed, upstream PR #21 updated → T03 done. geepers tiers approved → T12 done; geepers bug in `download_data_files` (kwargs into thread_map) fixed on geepers PR #3. S3 approved → first UNR snapshot uploaded (T48.2). VLM: products VLM-S1 / VLM-NI, one repo proposed, name pending. Gamma left: T04 docs, T05 URLs/version, T06 tag + delivery.
- 2026-10-07 T26 mostly done on stacked cal-disp PRs #3–#5: mask_file applied; reference-pixel rule ported so the core drops opera-utils[disp] (no dask/zarr, golden unchanged); dependency + peak-memory budget with CI step. geepers: grid downloads now atomic and verified (PR #4) so cal-disp can switch to it at the next pin bump.
- 2026-10-09 T37 wired: cal-disp calls `calibrate_pair`; schema-v2 options exposed; gamma golden passes at 1e-6; v0.5 on the golden pair closes exactly, 55 s / 7.27 GB. Wiring exposed two Venti memory regressions, fixed (`f7da024`, `2b2678b`, −3 GB). Next: T37.3 metadata + golden update, then T38 e2e on the four frames.
- 2026-10-09 T21–T25 done in disp2vlm_validation (PR #1; T25.3 batch deferred to T46): pipeline, classes, gate, reports, traceability, `run`/`summary` CLI; reproduces the trade-study F08882 numbers to the printed digit. T38 done: v0.5 (cal-disp PR #6) PASSES on all four frames vs gamma (sill −65 to −82%); Houston bias +1.14 mm/yr misses the 1 mm/yr target; bias issue #2 (inside the products, not chaining). Next: T39 (TS-G1), T42 (TS-U1: per-region offsets), T44 (TS-B1), T37.3 with T49.
- 2026-10-09 (later) T39 TS-G1, T43 TS-S1, T44.2 TS-B1 done (disp2vlm_validation `studies/TS-*/REPORT.md`). k per frame in the frame table (Venti PR #20). Open owner decisions: an additive σ₀ (2.6–4.9 mm) in σ_CAL; the class rule for 'coastal'. Issue #2 corrected: coherence power does not move the bias; the 50 km surface does. Next: T40 (stage frames 5–8), T42 (TS-U1), T47 (memory).
- 2026-10-09 (later) T42 TS-U1 phase 1 done: no-go for v0.5 (flag stays off). The v0.5 unwrap path had never run: fixed in cal-disp PR #7 and venti-dev PR #21. Galveston's real error is sub-region (T62); LA islands are beyond the 12 km anchor. Gate note: |bias| rule penalises fixing negative outliers on a positively biased frame (owner decision).
- 2026-10-09 (later) T47: memray-guided cuts, products bit-identical: v0.5 peak 7.27 → 4.0–4.1 GB, gamma 6.0 → 3.6–4.0 GB; budget 7.0 → 4.5 GB; worker_settings cap BLAS threads; 2-core run 1 min 43 s → t3.large fits. Unwrap correction stays off until a trade study confirms it (owner; note in T62).
