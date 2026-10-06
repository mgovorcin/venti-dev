# Venti

**Integrate GNSS calibration with InSAR products for accurate vertical land motion estimation**

---

## Overview

Venti is a Python toolkit designed to fuse GNSS data with OPERA-DISP product to calibrate and project InSAR line-of-sight measurements into vertical displacement estimates.

---

## Features

- GNSS calibration of InSAR products
- Projection from line-of-sight (LOS) to vertical displacement
- Line-of-sight (LOS) Decomposition

---

## Installation

1. Get the code:
```bash
git clone https://github.com/opera-adt/Venti.git
cd Venti
```
2. Create the environment (Venti's dependencies plus test and notebook tools):
```bash
mamba env create --file environment.yml
conda activate venti
```
3. Install Venti in editable mode:
```bash
python -m pip install -e .
```

Or, with [pixi](https://pixi.sh) (the lock file is committed, so this is the
reproducible route used in CI and by the operational image):
```bash
pixi install -e dev      # Venti + test, lint and docs tools
pixi run -e dev test     # pytest
pixi run -e dev lint     # pre-commit: ruff, black, mypy, nbstripout, SPDX headers
pixi run -e docs docs    # mkdocs build --strict
```
The `ops` environment (`pixi install -e ops`) is the default dependency set
only, i.e. what a production image installs.

With plain pip instead of conda, `python -m pip install -e ".[test]"` installs
Venti with its test tools. GDAL must already be available with NumPy support
(`python -c "from osgeo import gdal_array"` must work); to build the bindings
against a system GDAL, run `pip install numpy setuptools wheel` and then
`pip install --no-build-isolation "gdal==$(gdal-config --version)"`. The `notebooks`
extra adds notebook tools, and `analysis` adds the packages used by the GNSS
validation scripts (`numba`, `geepers`), which Venti itself does not need.

---

## Data Staging

Before running any calibration or decomposition workflow, you need to download and
prepare all input data for the target frame. The staging workflow handles this end-to-end.

### What gets staged

For a given OPERA frame ID and secondary date, staging prepares:

| Step | Output | Location |
|------|--------|----------|
| 1. Download DISP-S1 product | OPERA NetCDF interferogram | `<output_dir>/disp_s1/` |
| 2. Generate DEM | GLO30 GeoTIFF in WGS84 and native UTM | `<output_dir>/dem/` |
| 3. Generate LOS geometry | LOS ENU raster + incidence angle | `<output_dir>/los/` |
| 4. Tropospheric corrections | HRRR-based correction GeoTIFFs | `<output_dir>/tropo/` |
| 5. UNR GNSS velocities | Per-station `.tenv8` files + `velocities.parquet` | `<output_dir>/gnss/` |

DEM and LOS outputs are idempotent — re-running staging for a new date on the same
frame reuses existing files.

### Single product

Stage one DISP-S1 product by its secondary date.

**CLI:**

```bash
cd scripts/staging

# Stage all data for a single frame and date
python stage_frame_cli.py --frame-id 8887 --date 2016-06-15

# Custom output directory
python stage_frame_cli.py --frame-id 8887 --date 2016-06-15 \
    --output-dir /data/opera/frame_8887

# Skip optional steps
python stage_frame_cli.py --frame-id 8887 --date 2016-06-15 \
    --skip-tropo --skip-gnss

# Full options
python stage_frame_cli.py --help
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--frame-id` | required | OPERA frame identifier |
| `--date` | required | Secondary date (YYYY-MM-DD or YYYYMMDD) |
| `--output-dir` | `./staging` | Root directory for all outputs |
| `--num-workers` | `4` | Parallel workers for downloads |
| `--dem-buffer` | `10000.0` | Buffer in meters around frame for DEM |
| `--skip-tropo` | `False` | Skip tropospheric corrections |
| `--skip-gnss` | `False` | Skip GNSS download and velocity estimation |
| `--gnss-reference-frame` | `IGS20` | GNSS reference frame (`IGS20` or `IGS14`) |
| `--gnss-padding` | `0.0` | Extra meters beyond frame bounds for station search |
| `--gnss-start-year` | `2014.0` | Earliest year used for velocity estimation |

**Python API:**

```python
from pathlib import Path
from venti.workflow import run_data_staging

run_data_staging(
    frame_id=8887,
    date="2016-06-15",
    output_dir=Path("./data"),
)
```

### Time window (multiple products)

Stage all DISP-S1 products for a frame whose secondary date falls within a given
range.  DEM, LOS geometry, and GNSS velocities are produced once and shared
across all products.  Tropospheric corrections are batched over all unique epoch
sensing times.  The differential (secondary minus reference) is formed
per-product inside the calibration workflow.

**CLI:**

```bash
cd scripts/staging

# Stage all products for a frame in a given year
python stage_window_cli.py --frame-id 8887 --start 2016-01-01 --end 2016-12-31

# Custom output directory
python stage_window_cli.py --frame-id 8887 --start 2016-01-01 --end 2016-12-31 \
    --output-dir /data/opera/frame_8887

# Skip optional steps
python stage_window_cli.py --frame-id 8887 --start 2016-01-01 --end 2016-12-31 \
    --skip-tropo --skip-gnss

# Full options
python stage_window_cli.py --help
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--frame-id` | required | OPERA frame identifier |
| `--start` | required | Window start — secondary date (YYYY-MM-DD or YYYYMMDD) |
| `--end` | required | Window end — secondary date (YYYY-MM-DD or YYYYMMDD) |
| `--output-dir` | `./staging` | Root directory for all outputs |
| `--num-workers` | `4` | Parallel workers for downloads |
| `--dem-buffer` | `10000.0` | Buffer in meters around frame for DEM |
| `--skip-tropo` | `False` | Skip tropospheric corrections |
| `--skip-gnss` | `False` | Skip GNSS download and velocity estimation |
| `--gnss-reference-frame` | `IGS20` | GNSS reference frame (`IGS20` or `IGS14`) |
| `--gnss-padding` | `0.0` | Extra meters beyond frame bounds for station search |
| `--gnss-start-year` | `2014.0` | Earliest year used for velocity estimation |

**Python API:**

```python
from pathlib import Path
from venti.workflow import run_data_staging_window

disp_files = run_data_staging_window(
    frame_id=8887,
    start="2016-01-01",
    end="2016-12-31",
    output_dir=Path("./data"),
)
print(f"Staged {len(disp_files)} products")
```

### Preview available products before staging

Use `disp_cli.py` to check what products exist for a frame before downloading:

```bash
cd scripts/staging
python disp_cli.py preview --frame-id 8887 --start 2016-01-01 --end 2017-01-01 --print-dates
```

### Output directory structure

```
<output_dir>/
├── disp_s1/
│   ├── OPERA_L3_DISP-S1_IW_F08887_VV_20160101T*_20160115T*_*.nc
│   └── ...                              # one file per secondary date (window staging)
├── dem/
│   ├── dem_frame_8887.tif               # WGS84
│   └── dem_frame_8887_epsg32610.tif     # native UTM
├── los/
│   ├── los_enu_frame_8887.tif           # 3-band ENU raster
│   ├── incidence_angle_frame_8887.tif
│   ├── los_east.vrt
│   ├── los_north.vrt
│   └── los_up.vrt
├── tropo/
│   ├── tropo_urls.txt
│   ├── cropped_tropo/
│   ├── tropo_corrections/               # per-epoch absolute delays
│   └── tropo_corrections_{epsg}/        # reprojected per-epoch corrections
│       ├── tropo_{ref_timestamp}_{epsg}.tif
│       └── tropo_{sec_timestamp}_{epsg}.tif
└── gnss/
    ├── grid_latlon_lookup.txt
    ├── stations/
    │   └── *.tenv8
    └── velocities.parquet
```

---

## Running calibration from CLI

### 1. Generate configuration templates

```bash
python -m venti config --output-dir configs/
```

This writes two files to `configs/`:
- `runconfig.yaml` — run-specific settings (input/output paths, workflow type)
- `algorithm_parameters.yaml` — algorithm defaults that rarely need changing

### 2. Edit `runconfig.yaml`

Fill in the section that matches your workflow. For calibration:

```yaml
calibration_input_group:
  input_files: path/to/displacement/files   # directory of OPERA DISP NetCDF files
  los_file:    path/to/los_vectors.tif      # 3-band GeoTIFF (east, north, up)
  water_mask:  path/to/water_mask.tif       # GeoTIFF (1=land, 0=water)

product_path_group:
  product_path: output/

primary_executable:
  workflow_name: calibrate                  # or 'decompose'
```

### 3. Run the workflow

Two modes are available depending on whether you want to process a full directory of files or a single epoch.

**Batch run** — calibrates all `.nc` files found in `input_files`:

```bash
# Serial (default)
python -m venti run --config-file configs/runconfig.yaml

# Process up to 4 files concurrently (recommended range: 2–4)
python -m venti run --config-file configs/runconfig.yaml --n-workers 4
```

**Single-file run** — calibrates one specified displacement file:

```bash
python -m venti run-single --config-file configs/runconfig.yaml \
    --disp-file /path/to/epoch_001.nc
```

With optional tropospheric corrections (provide both the reference- and secondary-date files):

```bash
python -m venti run-single --config-file configs/runconfig.yaml \
    --disp-file /path/to/epoch_001.nc \
    --tropo-ref-file /path/to/tropo/tropo_{ref_timestamp}_{epsg}.tif \
    --tropo-sec-file /path/to/tropo/tropo_{sec_timestamp}_{epsg}.tif
```

`algorithm_parameters.yaml` is auto-discovered from the same directory as `runconfig.yaml`, so keep both files together. Increase verbosity with `--log-level DEBUG` on either command.

**Logging** — progress is printed to the console (timestamp, module, level). To also keep a log file, set `log_file: path/to/run.log` in `runconfig.yaml`; this works for `venti run`, `venti run-single`, `run_workflow` and `CalibrationWorkflow` in Python. `calibrate_timeseries`, `run_data_staging` and the staging CLIs take a `log_file` / `--log-file` argument. If your own script or notebook already configures logging (e.g. `logging.basicConfig`), Venti messages go through your setup instead of being printed twice. When calling lower-level functions directly, call `venti.configure_logging()` once to see their messages.

**Temporary files** — during `venti run`, `venti run-single` or a `CalibrationWorkflow` run, every temporary file (Venti's own and GDAL's) goes to `<scratch_path>/tmp`, so neither the working directory nor `/tmp` needs to be writable. Set `product_path_group.scratch_path` in `runconfig.yaml` to a writable directory; the run stops at the start with a clear error if it is not.

---

## Running from Python

**Batch run:**

```python
from venti.workflow.calibration import CalibrationWorkflow
from venti.workflow.config import load_config

config = load_config("configs/runconfig.yaml")
workflow = CalibrationWorkflow(config=config)

# Serial (default)
state = workflow.run()

# Parallel — process 4 files concurrently
state = workflow.run(n_workers=4)

print(f"Processed: {state.n_files_processed} / {state.n_files_total}")
for f in state.output_files:
    print(f"  {f.name}")
```

**Single-file run:**

```python
from pathlib import Path

from venti.workflow.calibration import CalibrationWorkflow
from venti.workflow.config import load_config

config = load_config("configs/runconfig.yaml")
workflow = CalibrationWorkflow(config=config)

# Without tropospheric correction
state = workflow.run_single(disp_file=Path("/path/to/epoch_001.nc"))

# With per-epoch tropospheric corrections
state = workflow.run_single(
    disp_file=Path("/path/to/epoch_001.nc"),
    tropo_ref_file=Path("/path/to/tropo/tropo_{ref_timestamp}_{epsg}.tif"),
    tropo_sec_file=Path("/path/to/tropo/tropo_{sec_timestamp}_{epsg}.tif"),
)

print(f"Output: {state.output_files[0]}")
```

`load_config` auto-discovers `algorithm_parameters.yaml` from the same directory as `runconfig.yaml`. You can also build the configuration fully in Python — see the [calibration workflow notebook](notebooks/calibration_workflow.ipynb) for a complete example.

---

### Setup for contributing


We use [pre-commit](https://pre-commit.com/) to automatically run linting, formatting, and [mypy type checking](https://www.mypy-lang.org/).
Additionally, we follow [`numpydoc` conventions for docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).
To install pre-commit locally, run:

```bash
pre-commit install
```
This adds a pre-commit hooks so that linting/formatting is done automatically. If code does not pass the checks, you will be prompted to fix it before committing.
Remember to re-add any files you want to commit which have been altered by `pre-commit`. You can do this by re-running `git add` on the files.

Since we use [black](https://black.readthedocs.io/en/stable/) for formatting and [flake8](https://flake8.pycqa.org/en/latest/) for linting, it can be helpful to install these plugins into your editor so that code gets formatted and linted as you save.

### Running the unit tests

After making functional changes and/or have added new tests, you should run pytest to check that everything is working as expected.

First, install the extra test dependencies:
```bash
python -m pip install --no-deps -e .[test]
```

Then run the tests:

```bash
pytest
```

---

## Using Venti's building blocks in your own workflow

Each step also works on its own with plain NumPy arrays; no Venti
configuration or OPERA product is needed.

```python
from venti import estimate_calibration_surface
from venti.gnss import project_to_los
from venti.workflow.config import CalibrationOptions

# GNSS LOS field (mm) on any grid: a NetCDF, a GeoTIFF or (x, y) coordinates.
gnss_los_mm = project_to_los(los_east, los_north, los_up, (x, y), station_gdf)

# Calibration surface (m) from arrays; corrections are removed before the fit.
result = estimate_calibration_surface(
    disp_m,
    gnss_los_mm / 1000,
    valid_mask,
    ref_point=(row, col),
    window_size=1000,                # pixels
    corrections=[tropo_m, solid_earth_tide_m],
    options=CalibrationOptions(),    # same options as algorithm_parameters.yaml
    downsample_factor=6,
)
calibrated = disp_m - result.surface
```

| Module | What it provides |
|---|---|
| `venti.surface` | `estimate_calibration_surface`: the full calibration step on arrays |
| `venti.gnss` | UNR download/readers, station rates or epoch displacements, LOS projection on any grid, cached LOS fields |
| `venti.spatial` | Windowed surface fit (`SpatialProcessor`), grid coordinates, interpolation, masked-region filling, residual-region detection, resampling |
| `venti.filtering` | Moving-window plane fit and low-pass filters |
| `venti.unwrap` | Unwrapping-error correction |
| `venti.io` | GeoTIFF/NetCDF readers and writers, OPERA correction layers |
