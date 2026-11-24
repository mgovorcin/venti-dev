# ADR-0022: A regenerated UNR rate grid, with its own sigma, may be the calibration input

- **Status:** accepted (owner, 2026-11-23; grid v33u_creep3_canada chosen 2026-11-24)
- **Date:** 2026-11-13; updated 2026-11-23 and 2026-11-24 (grid v33u_creep3_canada)
- **Source:** UNR grid regeneration study (`01_OPERA/CAL/unr_grid_regen`, README
  and LEDGER there); method in geepers (opera-adt/geepers#98, #99)
- **Supersedes, once accepted:** ADR-0005 in part (the published interpolated
  grid as the only input, and "trust its sigma"). ADR-0004 (constant grid
  operational) and ADR-0016 (frozen snapshot) stand.

## Decision

The time-constant GNSS grid that calibration reads may be one we rebuild from
UNR station positions, delivered as a frozen snapshot in the same layout as
the published one (R-G4), instead of the published UNR Version 0.3 grid.

The rebuilt grid:

1. takes MIDAS station rates through geepers, where MIDAS also avoids large
   uncatalogued steps;
2. keeps GPS Imaging's median filter to decide which stations are outliers;
3. leaves the stations of deforming areas out of the vertical, so the grid
   carries the regional field only (remove-restore on the GNSS side, at the
   stations and not at the nodes as R-G5 does today);
4. interpolates the vertical by local kriging with the spatial structure
   function, measured from the stations, as the variogram
   (`geepers.structure_kriging`), about the Caron (2018) GIA model as
   background;
5. interpolates the horizontals by least-squares collocation (LSC) about
   rigid plate motion (Pacific on the Pacific side, North America
   elsewhere): one global LSC in the plate-boundary belts (Southern
   California, Central Valley, Alaska, Hawaii), moving-window LSC elsewhere,
   blended over 100 km ("split hybrid");
6. removes models of fault creep from the station horizontals and restores
   them at the nodes (a fault creeping from the surface to depth Dc): the
   San Andreas creeping section through Parkfield (28.7 mm/yr, 15 km), the
   Hayward (7.2 mm/yr, 12 km) and the Calaveras (8.8 mm/yr, 15 km, tapered
   out at its north end), and adds 346 nodes on both sides of those traces
   (2, 6 and 12 km, every 6 km) so the sampled grid carries the steps. The
   Gulf of California stays cut by a breakline;
7. covers all of Canada (43,446 lattice nodes, the UNR layout extended
   north; 43,792 with the fault nodes, appended after the official lookup
   lines);
8. carries a sigma per component: the LSC or kriging sigma, plus for the
   vertical the GIA model spread far from stations, a term near the areas
   left out, and a term for the local misfit of the stations.

It becomes the operational input only after the e2e gate (plan T46) passes
with it and the owner signs off. Until then the published grid stays the
input, and this record stays proposed.

The chosen grid is snapshot `unr_grid_opera-v33u_creep3_canada_IGS20_20261123`,
in `s3://opera-adt/opera-ancillary/unr-grid/` next to the official snapshot
(uploaded and verified 2026-11-24; `CURRENT` unchanged). The earlier
candidates `unr_grid_opera-v33u_creep_canada_IGS20_20261120` (creeping section
only) and `unr_grid_opera-v33j_IGS20_20261112` are there too.
`cal-vlm-batch stage-cal --grid-snapshot` calibrates frames with a snapshot;
since cal-vlm-batch#9 it also reads an extended lookup, so the fault nodes are
staged where a frame reaches them.

## Context

ADR-0005 took the published grid as it is and trusted its sigma. Rebuilding it
from station positions showed where that does not hold:

- 184 of 365 nodes on the Baja peninsula carry North American motion; every
  station there moves with the Pacific plate.
- 522 nodes have zero sigma, and the sigma has no term for the distance to
  the stations: its median vertical value is 0.23 mm/yr within 15 km of
  stations and 0.46 at 150 to 250 km, where two stations that far apart
  differ by 0.72 mm/yr.
- The weighted median makes the surface step wherever it changes station:
  843 pairs of neighbouring nodes differ by over 2 mm/yr horizontally by more
  than three times their sigma.
- Subsiding basins reach the grid. In sparse networks the median cannot
  reject them, because the basin stations are the local majority (central
  Mexico).

Options considered, scored on stations held out in whole 75 km blocks
(7,820 stations, North America):

| | Published method rebuilt (GPS Imaging) | Kriging, this record |
|---|---|---|
| Typical error E/N/U (mm/yr) | 0.45 / 0.45 / 0.78 | 0.45 / 0.43 / 0.78 |
| Error over sigma E/N/U (1 is calibrated) | 1.52 / 1.44 / 1.05 | 0.93 / 0.87 / 1.04 |
| Sharp neighbour pairs not covered by 3 sigma | 843 | 158 |

Kriging does not improve accuracy where stations are dense. It gives a
surface without the median's steps and a sigma that is close to calibrated.
Removing steps from the positions before fitting rates was tried and made
rates worse. Keeping GPS Imaging and smoothing its grid afterwards was not
built: the smoothing and the sigma would be two unrelated additions.

Calibration test, ten benchmark frames, cal-disp 0.5.8 settings, only the grid
files swapped (2026-11-12):

- On stations outside every area the grid leaves out, RMSE against GNSS is
  lower with the rebuilt grid on nine frames and equal on the tenth; the
  double-difference sill changes by -6% to +4%.
- Inside the Houston bowl calibrated InSAR is about 2.2 mm/yr short of GNSS
  with the rebuilt grid (published: +0.1). The published grid closes that gap
  because its rim nodes carry the bowl; with the bowl stations kept in the
  rebuilt grid the gap is -0.45. The intent is that InSAR closes the gap on
  top of the gridded surface (owner, 2026-11-12). With cal-disp and the grid
  using the same areas the gap is 1.62 (from 1.87 on the same classes), and
  stable-station RMSE rises from 3.72 to 3.99: the rest of the gap is not a
  mismatch of areas and is not explained yet.

Update 2026-11-23, candidate v33u_creep_canada. Held out in 75 km blocks
(7,997 stations, Canada included), NMAD east/north/up 0.38 / 0.38 / 0.78 mm/yr
and error over sigma (|z| 68%) 0.90 / 0.88 / 0.91: LSC improves the
horizontals over kriging (0.45 / 0.43), and the extra vertical terms bring
Alaska's |z| from 1.35 to 1.06. Against independent solutions, reduced to the
regional field (median within 30 km, the vertical only beyond 25 km of the
areas left out), the horizontals differ from JPL and MEaSUREs by 0.32-0.37
mm/yr (published grid 0.39-0.45). The vertical difference (about 0.4) is the
level at which those centres differ from each other (0.33-0.36), not a grid
error. The horizontal error is 0.7-1.0 of the vertical, not the half that the
rule of thumb for positions gives, so it is measured, not imposed.

Calibration gate, ten frames, cal-disp 0.5.8 settings, the 137 deforming
areas masked in cal-disp and restored, the grid files swapped: stable-station
RMSE equal to or lower than with the published grid on nine frames; Hawaii
sill -7.6% (v33j: +14.9%). The two failures are not the grid's: the Central
Valley all-station sill (+10%, the same with the published grid and the 137
areas; +0.2% with the deforming-area stations out of the sill) and Houston's
subsiding-station condition (5.32 mm/yr, published grid with the same mask
5.37). Over the creeping section (F11116, 2016-2023) the creep model and fault
nodes lower gate RMSE from 4.21 (published) and 4.20 (without them) to 4.07,
and the 31 stations within 15 km of the trace from 3.36 to 2.96 mm/yr RMS. The
cubic RBF that Venti samples the grid with rings by about 4 mm/yr at the step,
but the calibration surface is long-wavelength: the calibrated velocity
changes smoothly (at most 3.5 mm/yr near the fault, 0.1 beyond the section
ends). A linear RBF or the creep model restored after sampling would follow
the step more closely; neither is needed for calibration, and both would
change values (Venti todo if ever wanted).

## Consequences

- R-G1 changes: the grid input is a frozen snapshot that may be regenerated
  by us from UNR positions; `snapshot.json` says which, and the snapshot id
  reaches the product provenance as it does today.
- R-G5 changes: stations of deforming areas are left out when the grid is
  built. cal-disp's node re-interpolation then has little left to do, and
  the areas used by the grid and by cal-disp have to be the same database.
- The deforming-area database grows. The study drafted 128 areas from GNSS
  alone; their extents need InSAR and owner review before any becomes
  curated (ADR-0009 stands).
- ADR-0017: the grid term of the uncertainty can use the grid's own sigma;
  the per-frame `k_grid`, fitted against the published sigma, has to be
  refitted or dropped.
- Validation scores stable and deforming stations apart, with the areas the
  grid leaves out as the classifier. Benchmark conditions on the deforming
  class were met with a grid that carried the bowls and need a decision.
- geepers stays the single source (ADR-0010): the method is
  `geepers.structure_kriging` and the recipe is
  `geepers/scripts/regenerate-unr-grid.py` (both on opera-adt/geepers main),
  which reproduces the study's grid exactly. The step detection behind the
  study's station rates and the deforming-area drafts are still study inputs.
- The method departs from the DISP-CAL ATBD section 2.5 (Hammond's GPS
  Imaging). The ATBD and the note to W. Hammond need the same account.
- Sparse arcs (the Aleutians) carry real block motion the grid follows only
  as far as their few stations allow; sigma there is probably too small.
- Creep is modelled on the San Andreas, Hayward and Calaveras only; other
  creeping faults (Superstition Hills, Imperial, Maacama) are smoothed like
  any other gradient, and the fault trace used looks about 1 km off near
  Parkfield. The fault models make the grid better near the faults (held
  out, Hayward stations within 5 km 3.3 -> 1.7 mm/yr, Parkfield 4.2 -> 2.7)
  but do not change calibration: on F11116 and F36542 the grid with all three
  models, the San Andreas alone and the creeping section alone calibrate
  alike (gate RMSE 4.07-4.10 and 2.66-2.67 mm/yr), because the calibration
  surface is long-wavelength.
- The snapshot has 43,792 nodes, more than the published 28,492. Readers that
  assume the official lookup have to read the lookup shipped with the
  snapshot instead: frame staging does since cal-vlm-batch#9; in geepers the
  snapshot is read with `CustomGridSource.from_snapshot` (opera-adt/geepers#100),
  since `UnrGridSource` stays UNR's product and knows only UNR's nodes.
- Open items are listed in the study's `LEDGER.md` (21 as of this date).
