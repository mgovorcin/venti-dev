# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tests for venti.sensor (plan T18, PRD R-X1)."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr

from venti.sensor import NISAR, S1, UnsupportedSensorError, get_sensor, sensor_for_file

S1_NAME = "OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_20251027T005420Z.nc"
STATIC_NAME = "OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif"


@pytest.fixture
def disp_file(tmp_path: Path) -> Path:
    """A small DISP-S1-shaped NetCDF: displacement, masks, /corrections, wavelength."""
    ny, nx = 8, 6
    rng = np.random.default_rng(0)
    ds = xr.Dataset(
        {
            "displacement": (
                ["y", "x"],
                rng.standard_normal((ny, nx)).astype(np.float32),
            ),
            "recommended_mask": (["y", "x"], np.ones((ny, nx), dtype=np.uint8)),
            "water_mask": (
                ["y", "x"],
                np.r_[np.zeros((2, nx)), np.ones((ny - 2, nx))].astype(np.uint8),
            ),
            "temporal_coherence": (
                ["y", "x"],
                rng.uniform(0.3, 0.9, (ny, nx)).astype(np.float32),
            ),
        },
        coords={"y": np.arange(ny, dtype=float), "x": np.arange(nx, dtype=float)},
    )
    path = tmp_path / S1_NAME
    ds.to_netcdf(path, engine="h5netcdf")
    corr = xr.Dataset(
        {"solid_earth_tide": (["y", "x"], np.full((ny, nx), 0.001, dtype=np.float32))}
    )
    corr.to_netcdf(path, group="corrections", mode="a", engine="h5netcdf")
    with h5py.File(path, "a") as f:
        f.create_group("identification").create_dataset(
            "radar_wavelength", data=np.float32(0.05546576)
        )
    return path


class TestConstants:
    def test_cycle_is_half_the_wavelength(self):
        assert S1.cycle_m == pytest.approx(0.02773, abs=1e-5)
        assert NISAR.cycle_m == pytest.approx(0.1207, abs=1e-3)

    def test_registry(self):
        assert get_sensor("s1") is S1
        with pytest.raises(UnsupportedSensorError, match="not implemented yet"):
            get_sensor("NISAR")
        with pytest.raises(UnsupportedSensorError, match="Unknown sensor"):
            get_sensor("ALOS")


class TestFilenames:
    def test_parse_product(self):
        pid = S1.parse_filename(Path("/data") / S1_NAME)
        assert pid.sensor == "S1"
        assert pid.frame_id == 8882
        assert pid.mode == "IW"
        assert pid.polarization == "VV"
        assert pid.version == "1.0"
        assert pid.reference_date == datetime(2022, 1, 11, 0, 26, 51)
        assert pid.secondary_date == datetime(2022, 7, 22, 0, 26, 57)
        assert pid.production_date.year == 2025
        assert pid.span_days == pytest.approx(192.0, abs=0.01)

    def test_parse_static(self):
        sid = S1.parse_static_filename(STATIC_NAME)
        assert (sid.frame_id, sid.platform, sid.layer, sid.version) == (
            8882,
            "S1A",
            "line_of_sight_enu",
            "1.0",
        )
        assert sid.reference_date == datetime(2014, 4, 3)

    def test_sensor_for_file(self):
        assert sensor_for_file(S1_NAME) is S1
        nisar = "OPERA_L3_DISP-NISAR_F00001_20250101T000000Z_20250113T000000Z_v0.1_20250201T000000Z.nc"
        with pytest.raises(UnsupportedSensorError, match="NISAR product"):
            sensor_for_file(nisar)
        with pytest.raises(UnsupportedSensorError, match="does not match"):
            sensor_for_file("random_file.nc")
        with pytest.raises(UnsupportedSensorError, match="not a S1"):
            S1.parse_filename("random_file.nc")
        assert not S1.matches("random_file.nc")


class TestReaders:
    def test_wavelength_from_product(self, disp_file):
        assert S1.read_wavelength(disp_file) == pytest.approx(0.05546576, rel=1e-6)

    def test_wavelength_fallback_and_mismatch(self, tmp_path, caplog):
        bare = tmp_path / "bare.nc"
        xr.Dataset({"displacement": (["y", "x"], np.zeros((2, 2)))}).to_netcdf(
            bare, engine="h5netcdf"
        )
        with caplog.at_level("WARNING"):
            assert S1.read_wavelength(bare) == S1.wavelength_m
        assert "no /identification/radar_wavelength" in caplog.text
        with h5py.File(bare, "a") as f:
            f.create_group("identification").create_dataset(
                "radar_wavelength", data=0.2413
            )
        with pytest.raises(UnsupportedSensorError, match="far from the S1 nominal"):
            S1.read_wavelength(bare)

    def test_displacement_and_masks(self, disp_file):
        disp = S1.open_displacement(disp_file)
        assert disp.shape == (8, 6)
        assert disp.dtype == np.float32
        water = S1.open_mask(disp_file, "water_mask")
        assert water.dtype == bool
        assert water.sum() == 6 * 6
        assert S1.open_mask(disp_file, "recommended_mask").all()
        with pytest.raises(ValueError, match="not a S1 mask layer"):
            S1.open_mask(disp_file, "displacement")

    def test_corrections(self, disp_file):
        assert S1.available_corrections(disp_file) == ("solid_earth_tide",)
        set_ = S1.open_correction(disp_file, "solid_earth_tide")
        assert set_ is not None
        assert set_.shape == (8, 6)
        assert np.allclose(set_, 0.001)
        assert S1.open_correction(disp_file, "ionospheric_delay") is None
        with pytest.raises(ValueError, match="not a S1 correction layer"):
            S1.open_correction(disp_file, "tropo")

    def test_no_corrections_group(self, tmp_path):
        bare = tmp_path / "bare.nc"
        xr.Dataset({"displacement": (["y", "x"], np.zeros((2, 2)))}).to_netcdf(
            bare, engine="h5netcdf"
        )
        assert S1.available_corrections(bare) == ()
