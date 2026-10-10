# Venti

**Fuse GNSS with OPERA DISP: calibrate InSAR line-of-sight displacement to an
absolute reference frame and project it to vertical land motion.**

Venti is the science library under two operational products:

| Product | What it makes | Repository |
|---|---|---|
| **DISP-CAL** (`cal-disp`) | one calibration surface per DISP-S1 granule; users apply `DISP − calibration` | [mgovorcin/cal-disp](https://github.com/mgovorcin/cal-disp) |
| **VLM** | per-pair vertical (and east) displacement from calibrated DISP plus GNSS | planned |

The design is written down in the [product requirements](specs.md); the work is
broken into dependency-ordered tasks in the [implementation plan](plan.md) and
tracked in the [task tracker](todo.md). Decisions that constrain the code are
[ADRs](decisions/README.md).

## What Venti does

- **GNSS → LOS.** Sample the UNR interpolated velocity grid (through
  [geepers](https://github.com/mgovorcin/geepers)), buffer it around the frame,
  exclude curated deforming areas, and project E/N/U to the radar line of sight.
- **Calibration surface.** Fit the long-wavelength (> 50 km) difference between
  DISP and GNSS with a robust, coherence-weighted local-linear kernel; InSAR
  keeps everything shorter. Tropospheric and solid-Earth-tide terms are removed
  and restored as explicit components.
- **Unwrap-cycle correction** (gated by trade study TS-U1): whole-cycle shifts
  for islands and cut-off peninsulas, decided against GNSS.
- **Uncertainty.** σ_CAL from the grid σ (inflated per frame), the fit, and
  the corrections.
- **Decomposition** (planned): ascending + descending → east and vertical,
  with north fixed from GNSS; projection to vertical for single geometries.

## Install

```bash
git clone https://github.com/mgovorcin/venti-dev.git Venti && cd Venti
pixi install -e dev          # or: pip install -e ".[test]"
pixi run -e dev test
```

GDAL must be importable with NumPy support (`from osgeo import gdal_array`);
see the README for building the bindings against a system GDAL.

## Use

Venti is driven through `cal-disp run` in operations. Directly:

```bash
venti config --output-dir run/      # writes runconfig.yaml + algorithm_parameters.yaml
venti run run/runconfig.yaml
```

The [API reference](api.md) documents the library entry points.
