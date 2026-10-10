# Venti, DISP-CAL and VLM: PRD

| | |
|---|---|
| **Status** | Draft v0.1, 2026-10-05 |
| **Owner** | Marin Govorcin |
| **Delivery sign-off (cal-disp)** | Talib (assumed; to confirm) |
| **Development repos** | `mgovorcin/venti-dev`, `mgovorcin/cal-disp`, `mgovorcin/geepers`, new `mgovorcin/<validation>`, new `mgovorcin/<vlm>` |
| **Upstream** | `opera-adt/*`; PRs only at milestone boundaries |

---

## 1. Project overview

### 1.1 Problem

OPERA DISP-S1 is *relative* displacement. Long-wavelength signal (plate motion, GIA, regional tectonics) is poorly constrained by InSAR and is contaminated by residual tropo, orbit and reference effects. Users who need absolute motion, or vertical land motion (VLM), must tie DISP to GNSS themselves.

### 1.2 Solution

Calibrate DISP against an interpolated GNSS grid:

- **GNSS constrains signal at wavelengths longer than 50 km.**
- **InSAR keeps its high resolution at shorter wavelengths.**
- Local fast deformation (subsiding basins, coseismic, volcanic) is handled by remove-restore, so InSAR resolves its full spatial extent.

| Deliverable | Kind | Role |
|---|---|---|
| **Venti** | Science library | Sensor-agnostic InSAR–GNSS fusion: GNSS→LOS sampling, calibration surface, unwrap correction, tropo modes, remove-restore, LOS decomposition and projection |
| **DISP-CAL (`cal-disp`)** | Operational SAS + product | One calibration product per DISP-S1 granule: a per-pair surface the user subtracts from DISP. The CLI and the runconfig are frozen. |
| **VLM** | Operational SAS + product (new repo) | Per-pair vertical (and east, where possible) displacement from calibrated DISP plus GNSS E/N |
| **Validation package** | VnV tool (new repo) | Independent e2e GNSS validation, release gates and VnV reports for DISP-CAL and VLM |

### 1.3 Current state (2026-10-05)

- **cal-disp `gamma-release`** (`VLM/DISP_CAL/gamma_release/cal-disp`, 52 commits past upstream `main`):
  - The only end-to-end path that works: a windowed bilinear plane against the UNR grid, optional tropo and SET, NetCDF output with `calibration` and `calibration_std`.
  - 478 tests and the golden F08882 (2022-01-11→07-22). Peak memory is 6.35 GB.
  - Venti is pinned at `2a7e61f`.
  - The checkout in `CAL/cal-disp` is a stale stub that writes zeros.
- **Venti** (`upstream/main` = `2a7e61f`, PR #18 "modular"; the exact commit cal-disp pins; fork `mgovorcin/venti-dev`, clone `00_tools/src/Venti`):
  - Working modules: `gnss/` (reference grid → LOS, UNR), `spatial/` (gap filling, interpolation, resample), `filtering/` (moving-window plane), `surface.py` (the calibration surface cal-disp calls), `unwrap/`, `models/` (ITRF PMM, GIA), `io/`, `workflow/`. 13 test files.
  - Verified bugs on `main`: the unwrap corrector receives the full radar wavelength λ from `surface.py` (`correct_region_offset(..., wavelength=wavelength_m)`), but one LOS cycle is λ/2; `models/load_itrf.py:110` imports `.plate_motion.euler_pole` from a module that is not a package.
  - Decomposition is a stub; there is no LOS→vertical step; PMM/GIA are not wired into the workflow.
  - The tropo research lives on the unmerged `models` branch (22 commits behind `main`). The checkouts in `01_OPERA/VLM/Venti` and `Venti_work/Venti` are stale snapshots.
- **Trade studies** (`VLM/DISP_CAL/trade studies/`) validated an improved algorithm on F08882, F08886, F16940 and F08622:
  - two-pass, then a local-linear surface with a 50 km cutoff, weighted by coh⁸;
  - F08882: DD sill 257→81 mm², velocity RMSE 5.07→3.28 mm/yr.
- **Docs:** no ATBD. The Product Spec doesn't match the delivered product (`gamma_release/DOCUMENT_FIXES.md`).

### 1.4 Goals

1. **G1.** Ship gamma 0.3 now: frozen algorithm, golden rebuilt in Docker.
2. **G2.** Move the validated trade-study science into Venti behind flags, then ship **cal-disp v0.5 (CalVal)** after it passes the 8-frame e2e gate and the edge cases.
3. **G3.** Make Venti sensor-agnostic, so DISP-NISAR (cal-disp v2) is a new `SensorSpec`, not a rewrite.
4. **G4.** Build VLM on Venti, with GNSS realization consistent with DISP-CAL.
5. **G5.** Keep the operational image small: lean cores and optional extras.
6. **G6.** Make every value change provable: e2e gate, deliberate golden regeneration, a VnV report.

### 1.5 Non-goals (for now)

- Changing the DISP-CAL CLI or runconfig schema.
- Calibrating seasonal, transient or coseismic signal with GNSS. The operational grid type is `constant` (secular).
- Unwrap correction inland. TS-U1 phase 2 comes later.
- Fixing any plate other than the per-frame lookup (no PA in California).
- Deriving remove-restore areas at run time. They are curated offline.

---

## 2. Core requirements

### 2.1 Product assumptions (must appear in the product spec)

| ID | Assumption |
|---|---|
| A1 | **Wavelength split.** The calibration contains GNSS-constrained signal at wavelengths longer than about 50 km. Shorter-wavelength InSAR signal is not modified, except by the explicit correction components (§3.1). |
| A2 | **Secular GNSS.** The operational constraint is the UNR **constant-velocity** grid × Δt. Seasonal, postseismic, coseismic and instrumental signal stays in calibrated DISP. Extrapolation past the end of the grid's data span is supported. |
| A3 | **Long-wavelength geophysics** (plate motion, GIA, regional tectonics) reaches the product through the GNSS grid. It is not modelled separately in the calibration. |
| A4 | **Grid resolution.** The UNR v0.3 grid uses GPS-Imaging median filtering with 25 km spacing, so it doesn't resolve wavelengths shorter than about 50 km. Grid σ is trusted as the weight, with the per-frame inflation k. |
| A5 | **Remove-restore.** Defo and event areas are excluded from the fit. InSAR keeps their full extent. |
| A6 | **Frame.** After calibration DISP is absolute in **IGS20 (ITRF2020)**, with no plate fixed. |

### 2.2 Interface requirements

| Surface | Rule |
|---|---|
| CLI (`cal-disp config/run/validate/validate-golden/download`, `opera_cal-disp`) | **Frozen** |
| Runconfig (`cal_disp_workflow`, `extra=forbid`) | **Frozen, byte-compatible with gamma** |
| `algorithm_parameters.yaml` | New fields only: additive, versioned, with defaults that reproduce gamma behaviour |
| Product | Additive only. Existing names, dtypes and filename are unchanged. New layers come through output packaging. |

- Remove-restore uses the existing `static_ancillary_group.defo_area_db_json` and `event_db_json`.

### 2.3 GNSS requirements

| ID | Requirement |
|---|---|
| R-G1 | Grid input is UNR v0.3, IGS20, through geepers. `constant` is the default. `variable` is experimental and allowed **only in reprocessing**, because it contains coseismic, postseismic, hydrological and instrumental signal. |
| R-G2 | GNSS buffer of 50 km outside the frame, with LOS extrapolated for nodes outside the swath |
| R-G3 | E/N/U→LOS projection with DISP-S1-STATIC LOS ENU |
| R-G4 | A **frozen grid snapshot** per release, mirrored to controlled storage. Its version and hash are recorded in the runconfig ancillary paths and the product metadata. Roll-forward is deliberate and gated by e2e. |
| R-G5 | GNSS-side remove-restore: drop grid nodes inside defo/event areas and re-interpolate them from surrounding nodes with GPS Imaging |
| R-G6 | Plate lookup per frame: **NA** default, **PA** Hawaii, **CA** (Caribbean) Puerto Rico. California stays NA. |

### 2.4 Algorithm requirements

All are behind `algorithm_parameters.yaml` flags. Gamma behaviour stays reproducible.

| ID | Requirement |
|---|---|
| R-S1 | **Two-pass:** (1) robust frame-wide tie; (2) unwrap correction on DISP − CAL₁, shifts applied to raw DISP (if enabled); (3) final surface |
| R-S2 | Local-linear kernel with a physical **half-response cutoff** `cutoff_wavelength_meters` (default 50 km), not a window size |
| R-S3 | **Fill gaps first.** The surface is continuous everywhere and never 0 on masked cells. |
| R-S4 | No global quantile outlier mask. Use local robust weights × **coherence^p** (p = 8; filled pixels weight 0.02). |
| R-S5 | Exclude InSAR pixels inside defo/event GeoJSON areas |
| R-T1 | Tropo mode `off \| stratified \| full` per frame, gated by DEM relief p5–p95: off below about 0.3 km, stratified at about 1.5 km or more, 0.3–1.5 km untested |
| R-T2 | SET from DISP `/corrections` |
| R-U1 | Unwrap correction **off** and `cal_unwrap_shift = 0` until TS-U1 phase 1 passes |
| R-X1 | `SensorSpec` carries λ, the product reader, filename parsing, available corrections and static-layer readers |

### 2.5 Uncertainty requirements

| ID | Requirement |
|---|---|
| R-E1 | σ_CAL² = (k·σ_grid)² + σ_fit² + σ_tropo² + σ_ref². **k is per frame**, from TS-G1, fitted on a station split held out from validation stations. σ_CAL grows inside interpolated defo/event areas. Unwrap shifts are treated as exact. σ_CAL is the uncertainty of the calibration surface (owner, 2026-10-09); per-pair non-secular station motion (TS-G1: 2.6–4.9 mm) is documented and available as the opt-in `uncertainty.sigma_nonsecular_meters`. |
| R-E2 | The product documents that σ_CAL covers only the calibration term. Users add DISP noise at wavelengths shorter than 50 km; the placeholder is about 10 mm until TS-S1 finishes. |
| R-E3 | VLM σ_U and σ_E are propagated through the WLS, including σ_CAL, σ_DISP, σ_N and the temporal-interpolation σ |
| R-E4 | σ realism (std(z) ∈ [0.8, 1.25], 95% coverage between 90% and 98%) is **reported** in v0.5, not gated |

### 2.6 Operational requirements

| ID | Requirement |
|---|---|
| R-O1 | Trigger: one DISP-CAL per DISP-S1 granule (forward processing), plus bulk reprocessing |
| R-O2 | Latency ≤ 72 h after all inputs are available |
| R-O3 | Runs on a **small EC2 instance**. No compute budget yet; the target is ≤ 4 GB peak, to be confirmed. |
| R-O4 | Deterministic: the same inputs, snapshot and config give the same product. Golden tolerance is 1e-6. |
| R-O5 | The golden is built **inside the Docker image** (tropo reprojection depends on the GDAL version) |
| R-O6 | The image excludes dask, zarr, matplotlib and jupyter. CI enforces an allow-list, a size budget and a peak-memory check. |

### 2.7 Validation requirements (release gate)

A run FAILS if any of these holds:

- the calibrated DD sill is not below raw;
- the sill rises more than 5% vs baseline and the rise is significant;
- velocity RMSE or |bias| rises more than 0.2 mm/yr;
- the golden-epoch check fails.

The release target is |bias| ≤ 1 mm/yr per frame. Validation uses **only independent GNSS** (UNR daily stations + MIDAS), never the calibration grid.

**CalVal v0.5 benchmark:**

| # | Category | Frame | Category-specific condition |
|---|---|---|---|
| 1 | Coastal subsidence + islands | F08882 Houston | subsiding stations: bias < 1, RMSE < 2 mm/yr |
| 2 | Stable inland | F08886 Oklahoma City | stable stations: \|bias\| < 1 mm/yr |
| 3 | High relief + islands | F16940 Los Angeles | sill drops; tropo mode doesn't make it worse |
| 4 | Northeast / GIA | F08622 New York | vertical pattern vs GNSS |
| 5 | Fast basin (defo area) | Central Valley (TBD) | continuous surface; stations outside unaffected; full extent of the bowl kept |
| 6 | Coseismic event | Ridgecrest 2019 (TBD) | no step leaks into CAL |
| 7 | Sparse GNSS | Great Basin or Montana (TBD) | realistic σ; no artifacts |
| 8 | Non-NA plate | Hawaii (PA) / Puerto Rico (CA) (TBD) | correct plate layer; Kīlauea as a defo area |

### 2.8 Edge cases (each needs a test or a benchmark frame)

| Edge case | Expected behaviour |
|---|---|
| Fast-subsiding basin | InSAR and GNSS nodes inside the area are excluded; the surface is interpolated through from stable ground |
| Coseismic event inside the pair | The step never leaks into CAL. Pairs that don't span the event behave as if there were no event DB. |
| Islands and cut-off peninsulas | Water-mask regions. Shifts are applied only after TS-U1. Free offsets can also remove real motion. |
| Plate boundary | NA-fixed DISP shows relative plate motion (about 40–50 mm/yr in California). This is expected. |
| Sparse GNSS | σ_CAL grows; no artifacts fitted on the interpolated grid |
| Masked or water cells | Surface is continuous, never 0 |
| Secondary date past the grid's data span | Constant-velocity extrapolation, recorded in metadata |
| Relief of 0.3–1.5 km | Tropo mode from the frame table; a risk until tested |

---

## 3. Core features

### 3.1 DISP-CAL product

Users apply it as **`DISP_cal = DISP − calibration`**.

| Layer | Content | Release |
|---|---|---|
| `calibration` | **Exact sum of all applied components** | gamma (frozen name and meaning) |
| `calibration_std` | σ_CAL (R-E1) | gamma (meaning upgraded in v0.5) |
| `cal_gnss_surface` | GNSS surface longer than 50 km, incl. reference-point tie | packaging |
| `cal_tropo` | tropo term, if applied | packaging |
| `cal_set` | solid Earth tide, if applied | packaging |
| `cal_unwrap_shift` | int cycles × λ/2; 0 until TS-U1 passes | packaging |
| `cal_reference_offset` | constant | packaging |
| `defo_area_mask` | remove-restore and event areas | packaging |
| `plate_motion` | rigid plate rotation (ITRF2020-PMM) → LOS × Δt; subtract it to get a plate-fixed product | packaging |
| `gnss_ve`, `gnss_vn` (+ σ) | GNSS **velocities** on the DISP grid, for VLM | packaging |

- **Invariant:** `calibration == Σ cal_*`, enforced by a test.
- `/metadata` records the flags for applied components, the algorithm parameters and the GNSS provenance (§4.3).

### 3.2 VLM product

| Item | Behaviour |
|---|---|
| Output | **Per-pair** vertical (and east, where asc+desc exist) displacement, with σ and a per-pixel mode flag |
| Asc+desc overlap | WLS for [E, U], with **N fixed from the GNSS grid** |
| Single geometry | `U = (LOS − e·E_gnss − n·N_gnss) / u`, flagged as valid only for long-wavelength horizontal motion |
| Temporal alignment | Moving-window denoising and temporal interpolation resample asc and desc onto common epochs |
| GNSS E/N | **Both options:** DISP-CAL `gnss_ve/vn` when present, otherwise the Venti `sample_gnss_enu`. The provenance hash must match, or the run fails. |
| Context | GIA (ICE-6G_D, Caron 2018) and PMM as context layers, not as corrections |

### 3.3 Calibration engine (Venti)

- Two-pass surface with a physical cutoff, gap filling and coherence-weighted robust fitting.
- GNSS buffer, LOS extrapolation and grid-node exclusion with GPS Imaging re-interpolation.
- Per-frame tropo modes and SET.
- Unwrap correction (gated): water-mask watershed regions, jumps measured against the anchored neighbour across water, a GNSS-direction veto, whole λ/2 cycles, and a per-region decisions CSV.
- σ_CAL with per-frame k.
- Sensor abstraction (S1 now, NISAR next).

### 3.4 Validation and VnV

- Metrics:
  - DD semivariogram sill (pairs ≥ 50 km, station-bootstrap CI);
  - velocity RMSE, bias and correlation vs MIDAS per station class (stable / subsiding / coastal / defo area / near event);
  - golden-epoch check;
  - VLM U vs GNSS Up and E vs GNSS East;
  - σ realism.
- Modes: candidate vs raw; candidate vs baseline; across frames.
- Outputs:
  - `metrics.json` and `stations.csv`;
  - an HTML/PDF VnV report per frame plus a release summary;
  - a requirements traceability table.
  - Audience: internal team plus the OPERA CalVal review.

---

## 4. Core components

### 4.1 Repository map and dependency chain

```
geepers                       opera-utils          Venti                            cal-disp / VLM             validation pkg
core:      UNR grid/stations  [disp] [tropo]       core: io, geometry, GNSS→LOS,    frozen CLI/runconfig,      e2e DD sill, MIDAS,
[grid]:    GPS Imaging, Euler                      config, SensorSpec               product I/O, packaging,    VnV reports,
[analysis]: MIDAS, stats,                          [calibration] [decomposition]    staging, golden            uses geepers[analysis]
            dask, zarr, pandera                    [models] [research] [all]
[all]
```

### 4.2 Ownership of shared logic

| Function | Single home | Action |
|---|---|---|
| UNR grid and station access | geepers core | remove the duplicates in Venti `gnss/unr.py` and cal-disp `_stage_unr.py`/`_unr.py` |
| Plate motion (Euler, ITRF2020-PMM) | geepers `[grid]` | delete Venti `models/plate_motion.py` |
| GPS Imaging re-interpolation | geepers `[grid]` | add an exclusion-area option |
| Tropo | cal-disp, copied into Venti | parity test keeps them identical; data via opera-utils `[tropo]` |
| Calibration, unwrap, remove-restore, decomposition | Venti | port from the trade studies |
| e2e validation | validation package | build from `validate_stack.py` |
| GIA / PMM rate models | Venti `[models]` | context layers for VLM |

### 4.3 Key interfaces

- **`sample_gnss_enu(grid_cfg, defo_db, geometry) → (vE, vN, vU, σ)`** is the only code path for sampling the grid. DISP-CAL and VLM both call it.
- **GNSS provenance in DISP-CAL `/metadata`:** grid version, type and snapshot hash; defo/event DB versions; the buffer; a hash of the sampled field.
- **`SensorSpec`:** λ, reader, filename parser, correction layers, static layers.
- **Frame-parameter table** (versioned): frame → plate, tropo mode, k, benchmark category.
- **Defo/event DB:** curated, versioned GeoJSON shipped as ancillary files.

### 4.4 Packaging rules

- **At most about 4 extras per package**, each tied to a real consumer, plus `[all]`.
- **Operational image:** `venti[calibration]` + `geepers` (core + `[grid]`) + `opera-utils[tropo]`.

### 4.5 Reusable prototype code (`VLM/DISP_CAL/trade studies/`)

| Path | Contents |
|---|---|
| `surface_fitting/continuous_surface/improved/surface_v2.py` | `loclin_v2`, `loclin_w` |
| `surface_fitting/continuous_surface/gnss_extend.py` | buffer and LOS extrapolation |
| `surface_fitting/continuous_surface/defo_area.py` | remove-restore areas |
| `two_pass/estimator_c.py`, `unwrap_step.py`, `tropo_mode.py` | two-pass chain |
| `e2e_validation/validate_stack.py`, `calibrate_stack.py`, `stage_tropo.py` | e2e validation |
| `unwrap_error_correction/unwrap_prototype/harness.py`, `score.py`, `truth.csv` | unwrap estimator test bench |

---

## 5. App/user flow

### 5.1 Operational DISP-CAL run (one granule)

1. A new DISP-S1 granule arrives. PCM stages the inputs:
   - DISP NetCDF;
   - DISP-S1-STATIC LOS ENU and DEM;
   - the frozen UNR grid snapshot;
   - OPERA TROPO (if the frame's tropo mode is not off);
   - defo/event GeoJSON;
   - the frame-parameter table.
2. `cal-disp run runconfig.yaml` (frozen interface) loads `algorithm_parameters.yaml`.
3. Venti:
   - samples the GNSS grid (buffer, exclusion, GPS Imaging re-interpolation) and projects it to LOS;
   - loads DISP and masks it (NaN, recommended, water, defo/event);
   - removes tropo and SET.
4. Pass 1: the robust frame-wide tie gives CAL₁.
5. (When enabled after TS-U1) unwrap correction runs on DISP − CAL₁, and the shifts are applied.
6. Pass 2: local-linear surface with the 50 km cutoff and coh^p weights, gaps filled. This produces `cal_gnss_surface` and σ_CAL.
7. Assemble `calibration = Σ components` and write the NetCDF, the browse PNG, metadata and provenance.
8. Deliver within ≤ 72 h of the inputs becoming available.

### 5.2 End user

1. Download a DISP-S1 granule and the matching DISP-CAL product.
2. Compute `DISP_cal = DISP − calibration`, which is absolute in IGS20.
3. Optionally subtract `plate_motion` to get a plate-fixed product, and inspect or undo components (`cal_tropo`, …).
4. Use `calibration_std` plus DISP's own noise at wavelengths shorter than 50 km for error budgets.
5. Treat `defo_area_mask` regions as InSAR-resolved, with the surface interpolated through them.

### 5.3 VLM run

1. Inputs: calibrated asc and desc DISP stacks plus DISP-CAL products.
2. Moving-window denoising and temporal interpolation resample them onto common epochs.
3. GNSS E/N come from DISP-CAL layers or `sample_gnss_enu`, with the provenance check.
4. Per pixel:
   - asc+desc gives a WLS solve for [E, U] with N fixed;
   - a single geometry gives a projection to U.
5. Write per-pair U (E), σ, the mode flag and GIA/PMM context.

### 5.4 Developer flow (any value-changing change)

1. Open an issue (feature / bug / science) and create a branch on the fork.
2. Implement in Venti behind a flag, with tests in the same commit.
3. Run `pixi run test` and the golden. If values changed, run `pixi run e2e` on the benchmark frames through the validation package.
4. Open a small single-concern PR with the VnV summary attached. The owner reviews.
5. For a release, regenerate the golden on purpose in Docker, with a changelog. Bump the version; open the upstream PR to `opera-adt` at the milestone.

---

## 6. Techstack

| Area | Choice |
|---|---|
| Language | Python ≥ 3.11 |
| Numerics | numpy, scipy, numba, scikit-image (watershed) |
| Geo / raster | rasterio, rioxarray, xarray, pyproj, shapely, geopandas (kept out of the core where possible) |
| Config | pydantic (frozen runconfig with `extra=forbid`), YAML |
| CLI | tyro |
| Product I/O | NetCDF4/CF (h5netcdf), GeoTIFF |
| GNSS | geepers (UNR grid/stations, GPS Imaging, Euler, MIDAS) |
| OPERA data | opera-utils (`[disp]`, `[tropo]`, frame DB) |
| Environments | pixi (`default`, `dev`, `ops`); conda-lock for Docker |
| Containers | Docker BuildKit, multi-stage, non-root, reproducible |
| Compute | AWS EC2 (small instance) for operations; AWS Batch via batchkit for benchmark and VnV runs |
| Quality | pre-commit (ruff, mypy, nbstripout, SPDX/Apache-2.0 + NOTICE), pre-commit.ci, pytest |
| Docs | mkdocs (material, mkdocstrings, mkdocs-jupyter) |
| Process | GitHub issue/PR templates, conventional commits, the `00_tools` skills (`git-commit-pr`, `unit-tests`, `workflow-regression`, `docker-build`), a `CLAUDE.md` per repo, Claude Code as implementer with owner review |

---

## 7. Implementation plan

### 7.1 Milestones

| Milestone | Window | Content | Done when |
|---|---|---|---|
| **M0: Gamma 0.3** | week 1 | Config confirmed with Talib; golden rebuilt in Docker (needs a build host); URLs; document fixes; upstream PR; tag; `addArtFile.sh` package | tag + delivery package |
| **M1: Foundations** | weeks 1–2 | geepers core/`[grid]`/`[analysis]`/`[all]`; Venti lean core + extras, packaging fixes, `SensorSpec`, bug fixes (λ/2, broken imports, arguments); validation package skeleton; pre-commit/pixi/templates/`CLAUDE.md` in all repos | CI green everywhere |
| **M2: Science port** | weeks 2–3 | R-G2…R-G6, R-S1…R-S5, R-T1, R-E1, `sample_gnss_enu` + provenance in Venti behind flags; cal-disp wired to it | e2e gate passes on F08882, F08886, F16940, F08622 |
| **M3: Trade studies** | weeks 2–6 | TS-G1 first, TS-U1 phase 1, TS-S1; benchmark staging (4 new stacks, defo/event GeoJSONs) | §7.2 exit criteria |
| **M4: cal-disp v0.5 CalVal** | ~6–8 weeks | 8-frame gate + edge cases; Docker golden; VnV report; small-EC2 target | VnV report signed off |
| **M5: VLM v0.1** | ~3–4 weeks after M4 | new repo; decomposition/projection; temporal resampling; E/N with provenance | VLM e2e vs GNSS Up/East |
| **M6: cal-disp v2 NISAR** | when beta DISP-NISAR arrives | NISAR `SensorSpec`; TS-N1; NISAR benchmark | e2e gate on NISAR frames |
| Later | — | TS-U1 phase 2; output-packaging layers; σ realism as a gate; `variable` grid evaluation | — |

### 7.2 Trade studies

| ID | Question | Truth / method | Exit criterion |
|---|---|---|---|
| **TS-G1** Grid fidelity | How well does the interpolated grid match individual sites? What is k per frame? | Interpolate the grid onto the DISP grid (step 1 only) and compare with independent stations by station density; start from geepers `cross_validation` | per-frame k table; grid bias/RMSE map |
| **TS-S1** DISP noise | What is σ_DISP at wavelengths shorter than 50 km? | Structure function of DISP − GNSS; temporal coherence and inversion residuals as proxies | per-frame noise model replaces the 10 mm placeholder |
| **TS-U1** Unwrap, phase 1 | Does whole-cycle correction work? | **Islands and cut-off peninsulas with trusted GNSS**, water-mask regions | shift = GNSS offset on every truth island; no shifts on decoys; sill doesn't rise |
| **TS-U1** Unwrap, phase 2 | Can inland errors be detected? | DISP `timeseries_inversion_residuals` | defined after phase 1 |
| **TS-N1** NISAR iono | Does GNSS + our iono handling beat the delivered correction? | Undo and redo the correction with the phase-screen layer in the cube; calibrate with and without it | same e2e gate |
| TS-T1 Tropo | Mode at 0.3–1.5 km relief | extend `tropo_check` | refined frame table |
| TS-B1 Stable bias | +2.8 mm/yr: site effect or real? coh⁸ vs coh¹⁶ vs a threshold? | extend `stable_bias` | weighting confirmed |

### 7.3 Release and golden policy

- **Gamma 0.3 line:** algorithm and golden frozen; fixes only.
- **v0.5 CalVal:** one deliberate golden regeneration in Docker, with a changelog explaining every value difference vs gamma.
- **After v0.5:** every value change needs the e2e gate + a new golden + a minor version bump.
- All development on `mgovorcin/*`. Upstream PRs to `opera-adt/*` happen at milestones.

### 7.4 Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Unwrap estimator validated on only 14 cases | a wrong shift costs 2.8 cm | off until TS-U1; whole cycles only; GNSS veto |
| Free offsets remove real island motion or tropo across water | biased islands | not in v0.5; part of TS-U1 |
| Grid can't resolve 50 km where GNSS is sparse | artifacts | TS-G1; per-frame k; σ_CAL grows |
| Velocity bias vs MIDAS (+3.7 mm/yr; input DISP / DOLPHIN chaining?) | bias gate fails | investigate in M2; report upstream if confirmed |
| Stable-site bias (+2.8 mm/yr) | coherence weighting may be masking real signal | TS-B1 |
| Tropo at 0.3–1.5 km relief untested | tropo may make things worse | TS-T1; default to off in that range |
| No Docker build host | M0 and M4 blocked | secure a build host in week 1 |
| NISAR ionosphere in the band longer than 50 km | GNSS surface absorbs iono | TS-N1 |
| 2–3 week window | over-commitment | window covers only M0–M2 + TS-G1; CalVal at 6–8 weeks |

### 7.5 Open questions

1. Owners and sign-off for geepers changes, the VLM gate and the UNR grid snapshot (assumed: Marin; Talib for cal-disp delivery).
2. Benchmark frame IDs for categories 5–8.
3. Compute budget and instance type; final memory target.
4. Names for the validation and VLM repos.
5. Product version string (v0.1 vs v1.0) and `reference_document` (unresolved from the gamma documents).
6. Snapshot mirror location (S3) and roll-forward cadence (proposal: every 6 months).
7. Does the OPERA CalVal plan require specific sites beyond the 8-frame set?

---

## Appendix: decision log (PRD interview, 2026-10-05)

| # | Decision |
|---|---|
| D1 | Venti = general science library; DISP-CAL = operational product. |
| D2 | DISP-CAL delivers a per-pair surface the user subtracts; VLM is a separate operational repo. |
| D3 | IGS20 calibration frame; plate lookup NA, PA for Hawaii, CA for Puerto Rico, no PA in California. |
| D4 | `constant` grid operational; `variable` experimental, reprocessing only. |
| D5 | The UNR interpolated grid is the calibration input; trust its σ; TS-G1 checks it against sites; no separate support layer. |
| D6 | `calibration` = Σ components, with component layers. |
| D7 | Unwrap correction gated on TS-U1 (islands with GNSS first, inversion residuals inland later). |
| D8 | CLI and runconfig frozen; new layers through output packaging. |
| D9 | Defo and event areas: curated static GeoJSON; InSAR and GNSS nodes excluded; GPS Imaging re-interpolation planned. |
| D10 | geepers is the single source for UNR access, plate motion (Euler) and GPS Imaging; tropo stays in cal-disp and is copied into Venti; separate validation package. |
| D11 | Extras-based packaging with `[all]`; avoid over-segmenting. |
| D12 | Sensor-agnostic Venti; S1 first; NISAR = cal-disp v2 with TS-N1. |
| D13 | Gamma frozen; v0.5 = CalVal release adopting the trade-study algorithm behind flags. |
| D14 | 8-frame CalVal benchmark. |
| D15 | VLM: per-pair; asc+desc WLS with N from GNSS, projection otherwise; moving-window temporal resampling; E/N from DISP-CAL or the grid, with a provenance check. **Owner, 2026-10-09:** the official product stays per pair; an option produces a velocity-domain product as well (calibrate the chained velocity once: the chained per-pair calibration spreads station residuals, NMAD 1.5 → 2.9 mm/yr on F08882, velocity-domain 2.2; disp2vlm_validation issue #3). |
| D16 | One run per granule; ≤ 72 h latency; small-EC2 goal; frozen grid snapshot. |
| D17 | k per frame; σ_DISP from TS-S1; σ realism reported, not gated. |
| D18 | Engineering standards; Claude Code in the workflow; all development on the fork. |
| D19 | PRD lives as Markdown in Venti `docs/specs.md`. |
