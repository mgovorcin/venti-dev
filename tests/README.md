# Venti tests

All tests use synthetic data: no network access, no OPERA downloads.

## Running

```bash
conda activate venti          # sets GDAL_DATA/PROJ_DATA; calling the env's
                              # python directly without activating breaks GDAL
pytest                        # tests + docstring examples in src/ (pyproject addopts)
pytest tests/test_surface.py  # one file
pytest --no-cov -q            # faster, no coverage report
```

Warnings are turned into errors (`filterwarnings = ["error", ...]` in
`pyproject.toml`), so a new warning fails the suite on purpose.

## What each file covers

| File | Covers |
|---|---|
| `test_surface.py` | `estimate_calibration_surface`: corrections removed before the fit and added back, downsampling, event masks and buffers, automatic outlier/region detection, smoothing options, input checks |
| `test_workflow_run.py` | `CalibrationWorkflow.run` / `run_single` end to end on small OPERA-format NetCDFs (fake GNSS source): outputs, naming, calibration quality, 1 vs 2 workers, uncertainty weighting, georeferenced average coherence |
| `test_gnss.py` | UNR readers, `sec - ref` sign, epoch tolerance, per-grid-type station cache, LOS projection, `constant` vs `variable` giving the same LOS field |
| `test_spatial.py` | Plane fitting, windowed fit, interpolation, `fill_gaps` (in memory, real zeros kept, scoped GDAL exceptions, read-only working directory) |
| `test_unwrap_corrections.py` | Unwrapping-error correction, including island correction and pixels outside corrected regions |
| `test_workflow.py` | Configuration models, YAML round trips, templates, date utilities, resampling |
| `test_io.py`, `test_product.py` | GeoTIFF/NetCDF I/O and product writers |
| `test_staging.py` | Data staging (DISP, DEM, LOS, tropo, GNSS) with mocked downloads |
| `test_cli.py` | The `venti` CLI run as a real subprocess |
| `test_basic.py` | Import smoke tests |

## Conventions

- Regression tests name the bug they guard in their docstring; each was
  checked to fail when the bug is put back.
- Tests that need GDAL, rasterio or geopandas use `pytest.importorskip`.
- Tests that change process-wide state (logging handlers, `tempfile.tempdir`,
  environment variables) restore it, or run in a subprocess.
