# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""End-to-end run of `CalibrationWorkflow` on small synthetic OPERA products.

Everything except the GNSS download is real: NetCDF reading, solid Earth
tide from ``/corrections``, reference-point auto-selection from the average
temporal coherence, GNSS LOS caching, the windowed fit, and output writing.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

# the stack drivers read DISP through opera_utils.disp (zarr): research tier
pytest.importorskip("zarr", reason="venti[research] not installed")

rasterio = pytest.importorskip("rasterio")
xr = pytest.importorskip("xarray")
from rasterio.transform import from_origin  # noqa: E402

from venti.workflow import calibration  # noqa: E402
from venti.workflow.config import (  # noqa: E402
    AlgorithmParameters,
    CalibrationInputGroup,
    CalibrationOptions,
    ProcessingOptions,
    ProductPathGroup,
    RunConfig,
    VentiConfig,
)

NY, NX, POSTING = 60, 80, 30.0
X0, Y0 = 500_000.0, 4_500_000.0
CRS = "EPSG:32618"
EPOCHS = [  # (reference, secondary) acquisition times in OPERA filenames
    ("20200101T000000Z", "20200301T000000Z"),
    ("20200101T000000Z", "20200601T000000Z"),
    ("20200101T000000Z", "20200901T000000Z"),
]
LOS = (-0.6, -0.1, 0.78)


def _gnss_velocity_mm_yr() -> np.ndarray:
    xx = np.arange(NX)[None, :].repeat(NY, 0)
    return 3.0 + 0.05 * xx


def _decimal_year(stamp: str) -> float:
    from venti.workflow.utils import get_file_dates

    return get_file_dates(f"x_{stamp}_{stamp}.nc")[0]


def _write_products(root):
    from pyproj import CRS as ProjCRS

    x = X0 + POSTING / 2 + POSTING * np.arange(NX)
    y = Y0 - POSTING / 2 - POSTING * np.arange(NY)
    yy, xx = np.mgrid[:NY, :NX]
    rng = np.random.default_rng(0)
    coherence = np.full((NY, NX), 0.5, "float32")
    coherence[25:35, 30:45] = 0.95  # best reference area
    disp_dir = root / "disp"
    disp_dir.mkdir()
    truth = {}
    for i, (ref, sec) in enumerate(EPOCHS):
        name = f"OPERA_L3_DISP-S1_IW_F08622_VV_{ref}_{sec}_v1.0_20260101T000000Z.nc"
        dt = _decimal_year(sec) - _decimal_year(ref)
        gnss_m = _gnss_velocity_mm_yr() * dt / 1000.0
        ramp = 0.003 * (i + 1) * yy / NY  # orbital-like error to remove
        tide = (0.002 * xx / NX - 0.001 * i).astype("float32")
        disp = gnss_m + ramp + tide + 2e-4 * rng.standard_normal((NY, NX))
        coords = {"y": y, "x": x}
        spatial_ref = xr.DataArray(0, attrs={"crs_wkt": ProjCRS(CRS).to_wkt()})
        xr.Dataset(
            {
                "displacement": (("y", "x"), disp.astype("float32")),
                "temporal_coherence": (("y", "x"), coherence),
                "water_mask": (("y", "x"), np.ones((NY, NX), "uint8")),
                "spatial_ref": spatial_ref,
            },
            coords,
        ).to_netcdf(disp_dir / name)
        xr.Dataset({"solid_earth_tide": (("y", "x"), tide)}, coords).to_netcdf(
            disp_dir / name, group="corrections", mode="a"
        )
        truth[name] = gnss_m

    transform = from_origin(X0, Y0, POSTING, POSTING)
    profile = {
        "driver": "GTiff",
        "height": NY,
        "width": NX,
        "crs": CRS,
        "transform": transform,
    }
    with rasterio.open(
        root / "los.tif", "w", count=3, dtype="float32", **profile
    ) as dst:
        for band, value in enumerate(LOS, start=1):
            dst.write(np.full((NY, NX), value, "float32"), band)
    with rasterio.open(
        root / "mask.tif", "w", count=1, dtype="uint8", **profile
    ) as dst:
        dst.write(np.ones((NY, NX), "uint8"), 1)
    return disp_dir, truth


class _FakeGNSS:
    """Stands in for a downloaded GNSSReference (constant-rate grid)."""

    grid_type = "constant"

    def compute_velocity_los(self, **_):
        return _gnss_velocity_mm_yr().astype("float32")

    def compute_velocity_los_std(self, **_):
        return np.full((NY, NX), 0.5, "float32")


@pytest.fixture
def products(tmp_path, monkeypatch):
    disp_dir, truth = _write_products(tmp_path)
    monkeypatch.setattr(calibration, "setup_gnss_reference", lambda **_: _FakeGNSS())
    return tmp_path, disp_dir, truth


def _config(root, disp_dir, out_name, **options):
    return VentiConfig(
        run_config=RunConfig(
            calibration_input_group=CalibrationInputGroup(
                input_files=disp_dir,
                los_file=root / "los.tif",
                water_mask=root / "mask.tif",
            ),
            product_path_group=ProductPathGroup(
                product_path=root / out_name, scratch_path=root / "scratch"
            ),
        ),
        algorithm_parameters=AlgorithmParameters(
            processing_options=ProcessingOptions(cal_downsample_factor=2),
            calibration_options=CalibrationOptions(
                window_size_meters=20 * POSTING,
                posting_meters=POSTING,
                unwrap_error_correction=False,
                **options,
            ),
        ),
    )


def _surfaces(state):
    import rasterio as rio

    out = {}
    for path in state.output_files:
        with rio.open(path) as src:
            out[path.name.split("_calibration_surface")[0] + ".nc"] = (
                src.read(1),
                src.descriptions[0],
                path.name,
            )
    return out


def test_run_calibrates_every_epoch(products):
    root, disp_dir, truth = products
    config = _config(root, disp_dir, "out")

    state = calibration.CalibrationWorkflow(config=config).run()

    assert (state.n_files_total, state.n_files_failed) == (3, 0)
    out_dir = root / "out"
    assert (out_dir / "average_temporal_coherence.tif").exists()
    assert (out_dir / "gnss_los_velocity.npy").exists()
    assert not (root / "scratch" / "tmp").exists()  # temp files cleaned up

    for name, (surface, description, filename) in _surfaces(state).items():
        assert filename.endswith(
            "_calibration_surface_constant_igs20_downsample2_set_nowrap.tif"
        )
        assert "IGS20 frame" in description
        assert "solid Earth tide" in description
        with xr.open_dataset(disp_dir / name) as ds:
            disp = ds["displacement"].values
        # The fake GNSS field is already in LOS, so calibration should recover it.
        residual = (disp - surface) - truth[name]
        assert np.nanstd(residual) < 5e-4
        assert abs(np.nanmean(residual)) < 5e-4


def test_parallel_workers_give_identical_surfaces(products):
    root, disp_dir, _ = products
    serial = _surfaces(
        calibration.CalibrationWorkflow(config=_config(root, disp_dir, "serial")).run()
    )
    parallel = _surfaces(
        calibration.CalibrationWorkflow(config=_config(root, disp_dir, "parallel")).run(
            n_workers=2
        )
    )
    assert serial.keys() == parallel.keys()
    for name in serial:
        np.testing.assert_array_equal(serial[name][0], parallel[name][0])


def test_run_single_matches_batch(products):
    root, disp_dir, _ = products
    batch = _surfaces(
        calibration.CalibrationWorkflow(config=_config(root, disp_dir, "batch")).run()
    )
    name = sorted(batch)[1]
    single_state = calibration.CalibrationWorkflow(
        config=_config(root, disp_dir, "single")
    ).run_single(disp_dir / name)
    single = _surfaces(single_state)

    np.testing.assert_array_equal(single[name][0], batch[name][0])


def test_uncertainty_weighting_runs_end_to_end(products):
    root, disp_dir, _ = products
    config = _config(root, disp_dir, "weighted", weight_fit_by_gnss_uncertainty=True)

    state = calibration.CalibrationWorkflow(config=config).run()

    assert state.n_files_failed == 0
    assert (root / "weighted" / "gnss_los_velocity_std.npy").exists()


def test_find_event_mask_file(tmp_path):
    disp = tmp_path / "epoch.nc"
    assert calibration.find_event_mask_file(None, disp) is None
    assert calibration.find_event_mask_file(tmp_path, disp) is None
    (tmp_path / "epoch_flood_mask.tif").touch()
    assert (
        calibration.find_event_mask_file(tmp_path, disp)
        == tmp_path / "epoch_flood_mask.tif"
    )


def test_setup_gnss_reference_needs_products(tmp_path):
    with pytest.raises(FileNotFoundError, match="No NetCDF files"):
        calibration.setup_gnss_reference(
            io_reader=SimpleNamespace(), input_files=tmp_path, product_path=tmp_path
        )


def test_average_coherence_is_georeferenced(products):
    root, disp_dir, _ = products
    from venti.workflow.utils import compute_average_temporal_coherence

    path = compute_average_temporal_coherence(
        sorted(disp_dir.glob("*.nc")), root / "coh"
    )

    with rasterio.open(path) as src, rasterio.open(root / "los.tif") as los:
        assert src.transform == los.transform
        assert src.crs == los.crs


def test_calibrate_timeseries_passes_options_through(products):
    from venti.workflow.run import calibrate_timeseries

    root, disp_dir, _ = products
    state = calibrate_timeseries(
        input_dir=disp_dir,
        los_file=root / "los.tif",
        mask_file=root / "mask.tif",
        output_dir=root / "functional",
        downsample_factor=2,
        window_size_meters=20 * POSTING,
        posting_meters=POSTING,
        unwrap_error_correction=False,
        apply_solid_earth_tide_correction=False,
        n_workers=2,
    )

    assert (state.n_files_total, state.n_files_failed) == (3, 0)
    for path in state.output_files:
        assert path.name.endswith(
            "_calibration_surface_constant_igs20_downsample2_nowrap.tif"
        )
