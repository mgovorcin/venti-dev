# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Unit tests for product classes (CalProduct and VlmProduct)."""

import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

# Optional imports
try:
    from rasterio.transform import Affine

    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    Affine = None

try:
    import h5netcdf  # type: ignore[unused-ignore]  # noqa: F401

    HAS_H5NETCDF = True
except ImportError:
    HAS_H5NETCDF = False

from venti.io.product import CalProduct, VlmProduct


@pytest.mark.skipif(not HAS_H5NETCDF, reason="h5netcdf not installed")
class TestCalProduct:
    """Test cases for CalProduct class."""

    def test_cal_product_basic_instantiation(self):
        """Test basic CalProduct instantiation."""
        shape = (100, 100)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_cal.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
            )

            assert product.output_path == output_path
            assert product.calibrated_displacement.shape == shape
            assert product.gnss_east.shape == shape
            assert product.gnss_north.shape == shape
            assert product.gnss_up.shape == shape
            assert product.displacement_std is None
            assert product.x_coords is not None
            assert product.y_coords is not None

    def test_cal_product_with_std_deviations(self):
        """Test CalProduct with standard deviations."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01
        disp_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001
        gnss_e_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.0001
        gnss_n_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.0001
        gnss_u_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.0001

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_cal_std.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
                displacement_std=disp_std,
                gnss_east_std=gnss_e_std,
                gnss_north_std=gnss_n_std,
                gnss_up_std=gnss_u_std,
            )

            assert product.displacement_std is not None
            assert product.gnss_east_std is not None
            assert product.gnss_north_std is not None
            assert product.gnss_up_std is not None

    @pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
    def test_cal_product_with_georeferencing(self):
        """Test CalProduct with transform and CRS."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01

        transform = Affine(30, 0, 500000, 0, -30, 4000000)
        crs = "EPSG:32610"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_cal_geo.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
                transform=transform,
                crs=crs,
            )

            assert product.transform == transform
            assert product.crs == crs
            # Coordinates should be generated from transform
            assert len(product.x_coords) == shape[1]
            assert len(product.y_coords) == shape[0]

    def test_cal_product_shape_validation(self):
        """Test that CalProduct validates array shapes."""
        shape = (50, 50)
        wrong_shape = (40, 40)
        disp = np.zeros(shape)
        gnss_e = np.zeros(shape)
        gnss_n = np.zeros(wrong_shape)  # Wrong shape!
        gnss_u = np.zeros(shape)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.nc"

            with pytest.raises(ValueError, match=r"shape.*does not match"):
                CalProduct(
                    output_path=output_path,
                    calibrated_displacement=disp,
                    gnss_east=gnss_e,
                    gnss_north=gnss_n,
                    gnss_up=gnss_u,
                )

    def test_cal_product_std_shape_validation(self):
        """Test that CalProduct validates std array shapes."""
        shape = (50, 50)
        wrong_shape = (40, 40)
        disp = np.zeros(shape)
        gnss_e = np.zeros(shape)
        gnss_n = np.zeros(shape)
        gnss_u = np.zeros(shape)
        disp_std = np.zeros(wrong_shape)  # Wrong shape!

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.nc"

            with pytest.raises(ValueError, match=r"shape.*does not match"):
                CalProduct(
                    output_path=output_path,
                    calibrated_displacement=disp,
                    gnss_east=gnss_e,
                    gnss_north=gnss_n,
                    gnss_up=gnss_u,
                    displacement_std=disp_std,
                )

    def test_cal_product_to_dataset(self):
        """Test converting CalProduct to xarray Dataset."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01
        disp_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
                displacement_std=disp_std,
            )

            ds = product.to_dataset()

            # Check dataset structure
            assert isinstance(ds, xr.Dataset)
            assert "calibrated_displacement" in ds.data_vars
            assert "gnss_east" in ds.data_vars
            assert "gnss_north" in ds.data_vars
            assert "gnss_up" in ds.data_vars
            assert "displacement_std" in ds.data_vars
            assert "x" in ds.coords
            assert "y" in ds.coords

            # Check attributes
            assert ds.attrs["product_type"] == "CAL"
            assert ds.attrs["title"] == "Venti Calibration Product"
            assert "history" in ds.attrs

            # Check data
            assert np.allclose(
                ds["calibrated_displacement"].values, disp, equal_nan=True
            )
            assert np.allclose(ds["gnss_east"].values, gnss_e, equal_nan=True)

    def test_cal_product_to_dataset_no_std(self):
        """Test Dataset creation without standard deviations."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
            )

            ds = product.to_dataset()

            # Std variables should not be present
            assert "displacement_std" not in ds.data_vars
            assert "gnss_east_std" not in ds.data_vars
            assert "gnss_north_std" not in ds.data_vars
            assert "gnss_up_std" not in ds.data_vars

    def test_cal_product_write(self):
        """Test writing CalProduct to NetCDF file."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "calibration.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
                reference_datetime=datetime(2020, 1, 1),
                secondary_datetime=datetime(2020, 1, 15),
            )

            product.write()

            # Verify file was created
            assert output_path.exists()

            # Verify file can be read
            ds = xr.open_dataset(output_path, engine="h5netcdf")
            assert "calibrated_displacement" in ds.data_vars
            assert ds.attrs["product_type"] == "CAL"
            assert "reference_datetime" in ds.attrs
            assert "secondary_datetime" in ds.attrs
            ds.close()

    def test_cal_product_write_with_compression(self):
        """Test writing with different compression levels."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "compressed.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
            )

            product.write(compression="zlib", complevel=9)

            assert output_path.exists()

            # Verify file is readable
            ds = xr.open_dataset(output_path, engine="h5netcdf")
            assert "calibrated_displacement" in ds.data_vars
            ds.close()

    def test_cal_product_with_metadata(self):
        """Test CalProduct with custom metadata."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01

        metadata = {
            "frame_id": 20697,
            "processing_version": "1.0.0",
            "gnss_stations_used": 15,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "metadata.nc"

            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
                metadata=metadata,
            )

            ds = product.to_dataset()

            assert ds.attrs["frame_id"] == 20697
            assert ds.attrs["processing_version"] == "1.0.0"
            assert ds.attrs["gnss_stations_used"] == 15


@pytest.mark.skipif(not HAS_H5NETCDF, reason="h5netcdf not installed")
class TestVlmProduct:
    """Test cases for VlmProduct class."""

    def test_vlm_product_basic_instantiation(self):
        """Test basic VlmProduct instantiation."""
        shape = (100, 100)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_vlm.nc"

            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
            )

            assert product.output_path == output_path
            assert product.east_displacement.shape == shape
            assert product.north_displacement.shape == shape
            assert product.up_displacement.shape == shape
            assert product.east_std is None
            assert product.x_coords is not None
            assert product.y_coords is not None

    def test_vlm_product_with_std_deviations(self):
        """Test VlmProduct with standard deviations."""
        shape = (50, 50)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01
        east_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001
        north_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001
        up_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_vlm_std.nc"

            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
                east_std=east_std,
                north_std=north_std,
                up_std=up_std,
            )

            assert product.east_std is not None
            assert product.north_std is not None
            assert product.up_std is not None

    @pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
    def test_vlm_product_with_georeferencing(self):
        """Test VlmProduct with transform and CRS."""
        shape = (50, 50)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01

        transform = Affine(120, 0, 500000, 0, -120, 4000000)
        crs = "EPSG:32610"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_vlm_geo.nc"

            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
                transform=transform,
                crs=crs,
            )

            assert product.transform == transform
            assert product.crs == crs

    def test_vlm_product_shape_validation(self):
        """Test that VlmProduct validates array shapes."""
        shape = (50, 50)
        wrong_shape = (40, 40)
        east = np.zeros(shape)
        north = np.zeros(wrong_shape)  # Wrong shape!
        up = np.zeros(shape)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.nc"

            with pytest.raises(ValueError, match=r"shape.*does not match"):
                VlmProduct(
                    output_path=output_path,
                    east_displacement=east,
                    north_displacement=north,
                    up_displacement=up,
                )

    def test_vlm_product_to_dataset(self):
        """Test converting VlmProduct to xarray Dataset."""
        shape = (50, 50)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01
        east_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.nc"

            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
                east_std=east_std,
            )

            ds = product.to_dataset()

            # Check dataset structure
            assert isinstance(ds, xr.Dataset)
            assert "east" in ds.data_vars
            assert "north" in ds.data_vars
            assert "up" in ds.data_vars
            assert "east_std" in ds.data_vars
            assert "x" in ds.coords
            assert "y" in ds.coords

            # Check attributes
            assert ds.attrs["product_type"] == "VLM"
            assert ds.attrs["title"] == "Venti VLM Product"
            assert "history" in ds.attrs

            # Check data
            assert np.allclose(ds["east"].values, east, equal_nan=True)
            assert np.allclose(ds["north"].values, north, equal_nan=True)
            assert np.allclose(ds["up"].values, up, equal_nan=True)

    def test_vlm_product_write(self):
        """Test writing VlmProduct to NetCDF file."""
        shape = (50, 50)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "decomposition.nc"

            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
                reference_datetime=datetime(2020, 1, 1),
                secondary_datetime=datetime(2020, 1, 15),
            )

            product.write()

            # Verify file was created
            assert output_path.exists()

            # Verify file can be read
            ds = xr.open_dataset(output_path, engine="h5netcdf")
            assert "east" in ds.data_vars
            assert "north" in ds.data_vars
            assert "up" in ds.data_vars
            assert ds.attrs["product_type"] == "VLM"
            assert "reference_datetime" in ds.attrs
            ds.close()

    def test_vlm_product_with_metadata(self):
        """Test VlmProduct with custom metadata."""
        shape = (50, 50)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01

        metadata = {
            "frame_id": 20697,
            "inversion_method": "weighted_least_squares",
            "min_geometries": 2,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "metadata.nc"

            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
                metadata=metadata,
            )

            ds = product.to_dataset()

            assert ds.attrs["frame_id"] == 20697
            assert ds.attrs["inversion_method"] == "weighted_least_squares"
            assert ds.attrs["min_geometries"] == 2


@pytest.mark.skipif(not HAS_H5NETCDF, reason="h5netcdf not installed")
class TestProductIntegration:
    """Integration tests for product classes."""

    def test_cal_product_roundtrip(self):
        """Test writing and reading CalProduct."""
        shape = (50, 50)
        disp = np.random.randn(*shape).astype(np.float32)
        gnss_e = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_n = np.random.randn(*shape).astype(np.float32) * 0.01
        gnss_u = np.random.randn(*shape).astype(np.float32) * 0.01
        disp_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "roundtrip.nc"

            # Create and write product
            product = CalProduct(
                output_path=output_path,
                calibrated_displacement=disp,
                gnss_east=gnss_e,
                gnss_north=gnss_n,
                gnss_up=gnss_u,
                displacement_std=disp_std,
                reference_datetime=datetime(2020, 1, 1),
                metadata={"test_key": "test_value"},
            )

            product.write()

            # Read back
            ds = xr.open_dataset(output_path, engine="h5netcdf")

            # Verify data
            assert np.allclose(
                ds["calibrated_displacement"].values, disp, equal_nan=True
            )
            assert np.allclose(ds["gnss_east"].values, gnss_e, equal_nan=True)
            assert np.allclose(ds["displacement_std"].values, disp_std, equal_nan=True)

            # Verify metadata
            assert ds.attrs["product_type"] == "CAL"
            assert ds.attrs["test_key"] == "test_value"
            assert "2020-01-01" in ds.attrs["reference_datetime"]

            ds.close()

    def test_vlm_product_roundtrip(self):
        """Test writing and reading VlmProduct."""
        shape = (50, 50)
        east = np.random.randn(*shape).astype(np.float32) * 0.01
        north = np.random.randn(*shape).astype(np.float32) * 0.01
        up = np.random.randn(*shape).astype(np.float32) * 0.01
        up_std = np.abs(np.random.randn(*shape).astype(np.float32)) * 0.001

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "roundtrip_vlm.nc"

            # Create and write product
            product = VlmProduct(
                output_path=output_path,
                east_displacement=east,
                north_displacement=north,
                up_displacement=up,
                up_std=up_std,
                secondary_datetime=datetime(2020, 1, 15),
                metadata={"method": "test"},
            )

            product.write()

            # Read back
            ds = xr.open_dataset(output_path, engine="h5netcdf")

            # Verify data
            assert np.allclose(ds["east"].values, east, equal_nan=True)
            assert np.allclose(ds["north"].values, north, equal_nan=True)
            assert np.allclose(ds["up"].values, up, equal_nan=True)
            assert np.allclose(ds["up_std"].values, up_std, equal_nan=True)

            # Verify metadata
            assert ds.attrs["product_type"] == "VLM"
            assert ds.attrs["method"] == "test"

            ds.close()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
