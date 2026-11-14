# ADR-0022: A regenerated UNR rate grid, with its own sigma, may be the calibration input

- **Status:** proposed
- **Date:** 2026-11-13
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
4. interpolates the remaining station rates by local kriging with the
   spatial structure function, measured from the stations, as the variogram
   (`geepers.structure_kriging`), with breaklines only where rates are truly
   discontinuous: the Gulf of California and the creeping section of the San
   Andreas;
5. carries the kriging sigma in the sigma columns.

It becomes the operational input only after the e2e gate (plan T46) passes
with it and the owner signs off. Until then the published grid stays the
input, and this record stays proposed.

Until then it is an alternative to try, not a replacement: the current
candidate, snapshot `unr_grid_opera-v33j_IGS20_20261112`, is in
`s3://opera-adt/opera-ancillary/unr-grid/` next to the official snapshot, with
`CURRENT` unchanged, and `cal-vlm-batch stage-cal --grid-snapshot` calibrates
frames with it (same nodes as the official staging, only the rates swapped).

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
- Open items are listed in the study's `LEDGER.md` (20 as of this date).
