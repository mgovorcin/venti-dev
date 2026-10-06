# Benchmarks

Numbers that gate algorithm decisions, with the script that reproduces them.
Value-changing results also go to the changelog (`workflow-regression`
policy in `CONTRIBUTING.md`).

## F08882 golden pair at real GNSS stations (plan T30.4)

**Question.** Does the v0.5 surface (local-linear kernel, 50 km half-response
cutoff, gap fill first, coherence and robust weights, two passes) beat the
gamma estimator (600 km windowed plane with the global residual mask) on the
pair the golden dataset is built from?

**Method.** `scripts/benchmark_golden_pair.py` calibrates the golden pair
F08882 2022-01-11 → 2022-07-22 with `venti.calibration.two_pass.calibrate_pair`
for each option set and compares `DISP − calibration` with real UNR stations
(IGS20 daily positions, not the grid the surface is fitted to) in LOS. InSAR
per station is the median of the valid pixels in a 7 × 7 window; the GNSS
reference is either the actual motion between the two dates (mean position
±3 d around each acquisition, "pair", 83 stations) or the MIDAS rate × 0.526 yr
("midas", 123 stations). The station table is
`tests/data/f08882_golden_pair_stations.csv`, built by trade study
`surface_fitting/gnss_pair_validation.py`. All variants use the same GNSS
LOS field (sampled at the fit-grid centres from the staged UNR nodes, no
buffer, as cal-disp does), the DISP solid-earth-tide correction, factor-6
downsampling and the DISP reference pixel; troposphere is off for every
variant because the HRES correction over-corrects this pair (trade study
`tropo_check`, frame table: F08882 `tropo.mode = off`).

```bash
NUMPY_MADVISE_HUGEPAGE=0 python scripts/benchmark_golden_pair.py \
    --golden-dir ../cal-disp/test_golden --out docs/benchmarks_f08882.md --json bench.json
```

**Reference points from the trade study** (same metric, cal-disp runs):
gamma golden config with tropo on 36.9 mm pair RMSE; gamma 600 km tropo off
21.4; gamma 11 km without the global mask 12.9; the trade-study `loclin_v2`
prototype 12.8–12.9.

### Results (2026-10-06, Venti stack through PR #18 + the two fixes below)

LOS, mm. Residual = InSAR − GNSS at the stations; "pair" n = 83, "midas"
n = 123. Grid vs stations (the UNR grid's own fidelity at the station pixels):
midas −0.1 ± 0.8 mm, pair −2.8 ± 7.8 mm, identical to the trade study. Runtime
is one variant on the full 7733 × 9464 frame (factor 6 fit grid), 32 cores;
peak RSS of the whole script 9.1 GB.

| Variant | pair bias | pair RMSE | pair NMAD | midas bias | midas RMSE | midas NMAD | fit std | s |
|---|---|---|---|---|---|---|---|---|
| raw DISP (relative to its reference point) | (−5.6) | (35.0) | 32.7 | (−3.2) | (33.0) | 31.5 | | 0 |
| **gamma 600 km** (golden config, tropo off) | +1.1 | 21.4 | 19.2 | +4.0 | 20.3 | 17.4 | | 23 |
| **loclin 50 km, two-pass, unweighted** | +3.3 | **13.6** | 12.0 | +5.6 | **12.3** | 8.8 | 12.0 | 61 |
| loclin 50 km, two-pass, coh⁸ | +3.5 | 14.3 | 11.4 | +5.6 | 13.2 | 9.0 | 10.2 | 59 |
| loclin 50 km, two-pass, coh⁸ + robust | +3.5 | 14.3 | 11.4 | +5.6 | 13.2 | 9.1 | 9.7 | 62 |
| loclin 50 km, two-pass, unweighted, defo excluded | +2.9 | 13.7 | 12.0 | +5.3 | 12.4 | 8.7 | 12.5 | 76 |
| loclin 50 km, two-pass, coh⁸ + robust, defo excluded | +2.9 | 14.6 | 10.6 | +5.2 | 13.4 | 9.0 | 11.6 | 81 |

Reading the table:

- **The gamma path reproduces cal-disp.** 21.4 / 20.3 mm with bias +1.1 / +4.0
  are the trade-study numbers for "ours, 600 km, tropo off" to the decimal, so
  `calibrate_pair`'s `windowed_plane` route is a faithful port (T30.4
  expectation met; the first run gave 24.3 mm because the GNSS sigma was
  being used as fit weights, which cal-disp only does with
  `weight_fit_by_gnss_uncertainty`; fixed).
- **The 50 km local-linear surface takes the station RMSE from 21.4 to
  13.6 mm (pair) and 20.3 to 12.3 mm (MIDAS)**, NMAD 19.2 → 12.0 and
  17.4 → 8.8. The trade-study prototype `loclin_v2` reached 12.8–12.9 mm
  with its own gap filling and GNSS field; the port is within 0.7 mm of it.
- **Coherence and robust weights do not help the station RMSE on this pair**
  (+0.7 mm); they lower the fit residual (12.0 → 9.7 mm) because they
  down-weight the rural pixels the stations do not sample. The trade studies
  chose coh⁸ on the stack-level double-difference sill (TS-B1, plan T44),
  which this single-pair station metric does not measure; the choice stays
  with T44.
- **Excluding the Houston subsidence bowls (7.7 % of the frame) moves the
  station bias by only −0.4 mm** and leaves the RMSE unchanged: the stations
  sit on the coherent urban area the 50 km kernel already follows. The
  remove-restore machinery works (5.6 M pixels excluded, surface interpolated
  through them); its value is the product inside the bowls, not these
  stations.
- **A +2 mm bias relative to gamma** (+3.3 vs +1.1 pair, +5.6 vs +4.0 MIDAS)
  is the one number to carry into the e2e gate (T38) and TS-G1 (T39): a
  50 km surface follows the frame's long-wavelength residual where a 600 km
  plane cannot, and whether that residual is DISP or GNSS is what the
  multi-frame bias check decides.
- **Runtime**: 60–80 s per loclin variant on the full frame, down from 290 s
  before the pass-1 tie moved to a coarse grid (`TIE_MAX_PX`); the gamma
  windowed fit takes 23 s.

Open items from this benchmark: T44 (weights on the sill), T38/T39 (bias),
T47 (peak memory: the full-resolution float32 copies dominate).
