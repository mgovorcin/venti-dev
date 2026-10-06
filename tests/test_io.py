# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Unit tests for I/O utilities (raster module)."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

# Optional imports
try:
    import rasterio as rio
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

from venti.io import raster


@pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
class TestReadGeoTIFF:
    """Test cases for read_geotiff function."""

    def test_read_geotiff_basic(self):
        """Test reading a basic GeoTIFF file."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            tif_path = tmpdir / "test.tif"

            # Create synthetic GeoTIFF
            data = np.random.rand(50, 50).astype(np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)
            crs = "EPSG:32611"

            with rio.open(
                tif_path,
                "w",
                driver="GTiff",
                height=50,
                width=50,
                count=1,
                dtype=data.dtype,
                crs=crs,
                transform=transform,
            ) as dst:
                dst.write(data, 1)

            # Read the file
            read_data, geo_info = raster.read_geotiff(tif_path)

            assert read_data.shape == (50, 50)
            assert np.allclose(read_data, data)
            assert "transform" in geo_info
            assert "crs" in geo_info
            assert "nodata" in geo_info
            assert "bounds" in geo_info
            assert geo_info["transform"] == transform
            assert geo_info["crs"] == rio.CRS.from_string(crs)

    def test_read_geotiff_multiband(self):
        """Test reading specific band from multiband GeoTIFF."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            tif_path = tmpdir / "multiband.tif"

            # Create multiband GeoTIFF
            data1 = np.ones((50, 50), dtype=np.float32)
            data2 = np.ones((50, 50), dtype=np.float32) * 2
            data3 = np.ones((50, 50), dtype=np.float32) * 3
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)

            with rio.open(
                tif_path,
                "w",
                driver="GTiff",
                height=50,
                width=50,
                count=3,
                dtype=np.float32,
                crs="EPSG:4326",
                transform=transform,
            ) as dst:
                dst.write(data1, 1)
                dst.write(data2, 2)
                dst.write(data3, 3)

            # Read band 2
            read_data, _ = raster.read_geotiff(tif_path, band=2)
            assert np.allclose(read_data, data2)

            # Read band 3
            read_data, _ = raster.read_geotiff(tif_path, band=3)
            assert np.allclose(read_data, data3)

    def test_read_geotiff_with_nodata(self):
        """Test reading GeoTIFF with nodata value."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            tif_path = tmpdir / "nodata.tif"

            # Create GeoTIFF with nodata
            data = np.random.rand(50, 50).astype(np.float32)
            data[0, 0] = -9999
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)

            with rio.open(
                tif_path,
                "w",
                driver="GTiff",
                height=50,
                width=50,
                count=1,
                dtype=data.dtype,
                crs="EPSG:32611",
                transform=transform,
                nodata=-9999,
            ) as dst:
                dst.write(data, 1)

            # Read the file
            read_data, geo_info = raster.read_geotiff(tif_path)

            assert geo_info["nodata"] == -9999
            assert read_data[0, 0] == -9999


@pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
class TestWriteGeoTIFF:
    """Test cases for write_geotiff function."""

    def test_write_geotiff_basic(self):
        """Test writing basic GeoTIFF with transform and CRS."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "output.tif"

            data = np.random.rand(50, 50).astype(np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)
            crs = "EPSG:32611"

            raster.write_geotiff(
                data, output_path, transform=transform, crs=crs, nodata=np.nan
            )

            assert output_path.exists()

            # Verify written data
            read_data, geo_info = raster.read_geotiff(output_path)
            assert read_data.shape == data.shape
            assert np.allclose(read_data, data, equal_nan=True)
            assert geo_info["transform"] == transform

    def test_write_geotiff_multiband(self):
        """Test writing multiband GeoTIFF."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "multiband.tif"

            # 3-band data
            data = np.random.rand(3, 50, 50).astype(np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)
            crs = "EPSG:32611"

            raster.write_geotiff(data, output_path, transform=transform, crs=crs)

            assert output_path.exists()

            # Verify band count
            with rio.open(output_path) as src:
                assert src.count == 3
                assert src.read().shape == data.shape

    def test_write_geotiff_with_descriptions(self):
        """Test writing GeoTIFF with band descriptions."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "described.tif"

            data = np.random.rand(2, 50, 50).astype(np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)
            descriptions = ["East velocity", "North velocity"]

            raster.write_geotiff(
                data,
                output_path,
                transform=transform,
                crs="EPSG:32611",
                descriptions=descriptions,
            )

            # Verify descriptions
            with rio.open(output_path) as src:
                assert src.descriptions[0] == "East velocity"
                assert src.descriptions[1] == "North velocity"

    def test_write_geotiff_with_reference_geotiff(self):
        """Test writing GeoTIFF using reference file (GeoTIFF)."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            ref_path = tmpdir / "reference.tif"
            output_path = tmpdir / "output.tif"

            # Create reference file
            ref_data = np.ones((50, 50), dtype=np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)
            crs = "EPSG:32611"

            with rio.open(
                ref_path,
                "w",
                driver="GTiff",
                height=50,
                width=50,
                count=1,
                dtype=ref_data.dtype,
                crs=crs,
                transform=transform,
                nodata=-9999,
            ) as dst:
                dst.write(ref_data, 1)

            # Write using reference
            data = np.random.rand(50, 50).astype(np.float32)
            raster.write_geotiff(data, output_path, reference_file=ref_path)

            # Verify georeferencing matches
            _, ref_geo_info = raster.read_geotiff(ref_path)
            _, out_geo_info = raster.read_geotiff(output_path)

            assert out_geo_info["transform"] == ref_geo_info["transform"]
            assert out_geo_info["crs"] == ref_geo_info["crs"]

    def test_write_geotiff_with_reference_netcdf(self):
        """Test writing GeoTIFF using reference file (NetCDF)."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            ref_path = tmpdir / "reference.nc"
            output_path = tmpdir / "output.tif"

            # Create reference NetCDF
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            disp_data = np.random.rand(len(y), len(x))

            ds = xr.Dataset(
                {"displacement": (["y", "x"], disp_data)}, coords={"x": x, "y": y}
            )
            ds.to_netcdf(ref_path)

            # Write using NetCDF reference
            data = np.random.rand(len(y), len(x)).astype(np.float32)
            raster.write_geotiff(data, output_path, reference_file=ref_path)

            assert output_path.exists()

    def test_write_geotiff_masked_array(self):
        """Test writing masked array as GeoTIFF."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "masked.tif"

            # Create masked array
            data = np.random.rand(50, 50).astype(np.float32)
            mask = np.zeros((50, 50), dtype=bool)
            mask[0:10, 0:10] = True
            masked_data = np.ma.masked_array(data, mask=mask)

            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)

            raster.write_geotiff(
                masked_data, output_path, transform=transform, crs="EPSG:32611"
            )

            # Verify masked values are replaced with nodata
            read_data, _ = raster.read_geotiff(output_path)
            assert np.isnan(read_data[0, 0])

    def test_write_geotiff_missing_georef(self):
        """Test error when georeferencing is missing."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "output.tif"

            data = np.random.rand(50, 50).astype(np.float32)

            with pytest.raises(
                ValueError, match=r"transform and crs.*must be provided"
            ):
                raster.write_geotiff(data, output_path)

    def test_write_geotiff_custom_dtype(self):
        """Test writing with custom data type."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "int16.tif"

            data = (np.random.rand(50, 50) * 1000).astype(np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)

            raster.write_geotiff(
                data,
                output_path,
                transform=transform,
                crs="EPSG:32611",
                dtype="int16",
                nodata=-9999,
            )

            with rio.open(output_path) as src:
                assert src.dtypes[0] == "int16"

    def test_write_geotiff_compression(self):
        """Test writing with compression."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            output_path = tmpdir / "compressed.tif"

            data = np.random.rand(50, 50).astype(np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)

            raster.write_geotiff(
                data,
                output_path,
                transform=transform,
                crs="EPSG:32611",
                compress="deflate",
            )

            with rio.open(output_path) as src:
                # Compression values are uppercase in rasterio
                assert src.compression.value.lower() == "deflate"


@pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
class TestReadNetCDF:
    """Test cases for read_netcdf function."""

    def test_read_netcdf_basic(self):
        """Test reading NetCDF with basic structure."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "test.nc"

            # Create synthetic data
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            disp_data = np.random.rand(len(y), len(x))
            mask_data = np.ones((len(y), len(x)), dtype=bool)

            ds = xr.Dataset(
                {
                    "displacement": (["y", "x"], disp_data),
                    "water_mask": (["y", "x"], mask_data),
                },
                coords={"x": x, "y": y},
            )

            ds.to_netcdf(nc_path)

            # Test reading
            disp, mask, geo_info = raster.read_netcdf(nc_path)

            assert disp.shape == disp_data.shape
            assert mask.shape == mask_data.shape
            assert np.allclose(disp, disp_data)
            assert "transform" in geo_info
            assert "crs" in geo_info

    def test_read_netcdf_no_mask(self):
        """Test reading NetCDF without mask variable."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "no_mask.nc"

            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            disp_data = np.random.rand(len(y), len(x))

            ds = xr.Dataset(
                {"displacement": (["y", "x"], disp_data)}, coords={"x": x, "y": y}
            )

            ds.to_netcdf(nc_path)

            # Read without mask
            disp, mask, _ = raster.read_netcdf(nc_path, mask_variable=None)

            assert disp.shape == disp_data.shape
            assert mask is None

    def test_read_netcdf_custom_variable(self):
        """Test reading NetCDF with custom variable name."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "custom.nc"

            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            velocity_data = np.random.rand(len(y), len(x))

            ds = xr.Dataset(
                {"velocity": (["y", "x"], velocity_data)}, coords={"x": x, "y": y}
            )

            ds.to_netcdf(nc_path)

            # Read custom variable
            data, _, _ = raster.read_netcdf(nc_path, variable="velocity")

            assert np.allclose(data, velocity_data)

    def test_read_netcdf_with_fillvalue(self):
        """Test reading NetCDF defaults to nan when no fillvalue specified."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "fillvalue.nc"

            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            disp_data = np.random.rand(len(y), len(x))

            ds = xr.Dataset(
                {"displacement": (["y", "x"], disp_data)}, coords={"x": x, "y": y}
            )

            ds.to_netcdf(nc_path)

            # Read and check nodata defaults to nan
            _, _, geo_info = raster.read_netcdf(nc_path)

            # When no _FillValue or missing_value, defaults to nan
            assert np.isnan(geo_info["nodata"])

    def test_read_netcdf_missing_variable(self):
        """Test error for missing variable."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "missing.nc"

            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)

            ds = xr.Dataset(
                {"other": (["y", "x"], np.ones((len(y), len(x))))},
                coords={"x": x, "y": y},
            )

            ds.to_netcdf(nc_path)

            with pytest.raises(ValueError, match=r"Displacement variable.*not found"):
                raster.read_netcdf(nc_path)

    def test_read_netcdf_wrong_dimensions(self):
        """Test error for wrong data dimensions."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "3d.nc"

            # Create 3D data
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            t = np.arange(0, 10)
            data_3d = np.random.rand(len(t), len(y), len(x))

            ds = xr.Dataset(
                {"displacement": (["t", "y", "x"], data_3d)},
                coords={"x": x, "y": y, "t": t},
            )

            ds.to_netcdf(nc_path)

            with pytest.raises(ValueError, match="Expected 2D data"):
                raster.read_netcdf(nc_path)


@pytest.mark.skipif(not HAS_H5NETCDF, reason="h5netcdf not installed")
def test_read_netcdf_correction(tmp_path):
    path = tmp_path / "product.nc"
    xr.Dataset({"displacement": (("y", "x"), np.zeros((3, 4)))}).to_netcdf(path)
    no_group = raster.read_netcdf_correction(path, "solid_earth_tide")
    tide = np.arange(12.0).reshape(3, 4)
    xr.Dataset({"solid_earth_tide": (("y", "x"), tide)}).to_netcdf(
        path, group="corrections", mode="a"
    )

    assert no_group is None
    assert raster.read_netcdf_correction(path, "ionospheric_delay") is None
    np.testing.assert_array_equal(
        raster.read_netcdf_correction(path, "solid_earth_tide"), tide
    )


class TestUpdateNetCDFVariable:
    """Test cases for update_netcdf_variable function."""

    def test_update_existing_variable(self):
        """Test updating an existing variable in NetCDF."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "update.nc"

            # Create NetCDF with initial data
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            initial_data = np.ones((len(y), len(x)))

            ds = xr.Dataset(
                {"displacement": (["y", "x"], initial_data)}, coords={"x": x, "y": y}
            )

            ds.to_netcdf(nc_path, engine="h5netcdf")

            # Update the variable
            new_data = np.ones((len(y), len(x))) * 2.0
            raster.update_netcdf_variable(nc_path, "displacement", new_data)

            # Verify update
            ds_updated = xr.open_dataset(nc_path, engine="h5netcdf")
            assert np.allclose(ds_updated["displacement"].values, new_data)
            ds_updated.close()

    def test_update_with_attributes(self):
        """Test updating variable with new attributes."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "attrs.nc"

            # Create NetCDF
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            data = np.ones((len(y), len(x)))

            ds = xr.Dataset(
                {"displacement": (["y", "x"], data)}, coords={"x": x, "y": y}
            )

            ds.to_netcdf(nc_path, engine="h5netcdf")

            # Update with attributes
            new_data = np.ones((len(y), len(x))) * 2.0
            attrs = {"units": "meters", "description": "Updated displacement"}
            raster.update_netcdf_variable(
                nc_path, "displacement", new_data, variable_attrs=attrs
            )

            # Verify attributes
            ds_updated = xr.open_dataset(nc_path, engine="h5netcdf")
            assert ds_updated["displacement"].attrs["units"] == "meters"
            assert (
                ds_updated["displacement"].attrs["description"]
                == "Updated displacement"
            )
            ds_updated.close()

    def test_update_missing_variable_create_false(self):
        """Test error when variable doesn't exist and create_if_missing=False."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "missing.nc"

            # Create NetCDF without target variable
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            data = np.ones((len(y), len(x)))

            ds = xr.Dataset({"other": (["y", "x"], data)}, coords={"x": x, "y": y})

            ds.to_netcdf(nc_path, engine="h5netcdf")

            # Try to update non-existent variable
            new_data = np.ones((len(y), len(x))) * 2.0

            with pytest.raises(ValueError, match=r"Variable.*not found"):
                raster.update_netcdf_variable(
                    nc_path, "displacement", new_data, create_if_missing=False
                )

    def test_update_missing_variable_create_true(self):
        """Test error when trying to create new variable (not implemented)."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "create.nc"

            # Create NetCDF
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            data = np.ones((len(y), len(x)))

            ds = xr.Dataset({"other": (["y", "x"], data)}, coords={"x": x, "y": y})

            ds.to_netcdf(nc_path, engine="h5netcdf")

            # Try to create new variable (should raise NotImplementedError)
            new_data = np.ones((len(y), len(x))) * 2.0

            with pytest.raises(NotImplementedError, match="not yet implemented"):
                raster.update_netcdf_variable(
                    nc_path, "new_variable", new_data, create_if_missing=True
                )


@pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
class TestGetBounds:
    """Test cases for get_bounds function."""

    def test_get_bounds_geotiff(self):
        """Test getting bounds from GeoTIFF."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            tif_path = tmpdir / "bounds.tif"

            # Create GeoTIFF with known bounds
            data = np.ones((50, 50), dtype=np.float32)
            transform = Affine(10.0, 0.0, 100.0, 0.0, -10.0, 500.0)

            with rio.open(
                tif_path,
                "w",
                driver="GTiff",
                height=50,
                width=50,
                count=1,
                dtype=data.dtype,
                crs="EPSG:32611",
                transform=transform,
            ) as dst:
                dst.write(data, 1)

            # Get bounds
            S, N, W, E = raster.get_bounds(tif_path)

            # Expected bounds from transform
            # W = 100, E = 100 + 50*10 = 600
            # N = 500, S = 500 - 50*10 = 0
            assert W == 100.0
            assert E == 600.0
            assert S == 0.0
            assert N == 500.0

    def test_get_bounds_netcdf(self):
        """Test getting bounds from NetCDF."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "bounds.nc"

            # Create NetCDF with known coordinates
            x = np.arange(100, 200, 10)
            y = np.arange(500, 600, 10)
            data = np.ones((len(y), len(x)))

            ds = xr.Dataset(
                {"displacement": (["y", "x"], data)}, coords={"x": x, "y": y}
            )

            ds.to_netcdf(nc_path)

            # Get bounds
            S, N, W, E = raster.get_bounds(nc_path)

            assert W == 100.0
            assert E == 190.0
            assert S == 500.0
            assert N == 590.0

    def test_get_bounds_netcdf_missing_coords(self):
        """Test error when NetCDF has no x/y coordinates."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nc_path = tmpdir / "no_coords.nc"

            # Create NetCDF without x/y coords
            data = np.ones((10, 10))
            ds = xr.Dataset({"displacement": (["row", "col"], data)})

            ds.to_netcdf(nc_path)

            with pytest.raises(ValueError, match="No x/y coordinates"):
                raster.get_bounds(nc_path)


@pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
class TestExtractNetCDFGeoref:
    """Test cases for _extract_netcdf_georef helper function."""

    def test_extract_with_spatial_ref(self):
        """Test georef extraction with spatial_ref variable."""
        x = np.arange(0, 100, 10)
        y = np.arange(0, 50, 10)
        data = np.ones((len(y), len(x)))

        ds = xr.Dataset(
            {"displacement": (["y", "x"], data), "spatial_ref": ([], 0)},
            coords={"x": x, "y": y},
        )

        # Add CRS and GeoTransform
        ds["spatial_ref"].attrs["crs_wkt"] = (
            'PROJCS["WGS 84 / UTM zone 11N",GEOGCS["WGS'
            ' 84",DATUM["WGS_1984",SPHEROID["WGS'
            ' 84",6378137,298.257223563]],PRIMEM["Greenwich",0],'
            'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
            'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-117],'
            'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
            'PARAMETER["false_northing",0],UNIT["metre",1]]'
        )
        ds["spatial_ref"].attrs["GeoTransform"] = "0.0 10.0 0.0 50.0 0.0 -10.0"

        geo_info = raster._extract_netcdf_georef(ds, data.shape)

        assert "crs" in geo_info
        assert "transform" in geo_info
        assert geo_info["transform"] is not None

    def test_extract_from_coords(self):
        """Test georef extraction from x/y coordinates."""
        x = np.arange(0, 100, 10)
        y = np.arange(0, 50, 10)
        data = np.ones((len(y), len(x)))

        ds = xr.Dataset({"displacement": (["y", "x"], data)}, coords={"x": x, "y": y})

        geo_info = raster._extract_netcdf_georef(ds, data.shape)

        assert "transform" in geo_info
        assert geo_info["transform"] is not None

    def test_extract_default_crs(self):
        """Test default CRS when not specified."""
        x = np.arange(0, 100, 10)
        y = np.arange(0, 50, 10)
        data = np.ones((len(y), len(x)))

        ds = xr.Dataset({"displacement": (["y", "x"], data)}, coords={"x": x, "y": y})

        geo_info = raster._extract_netcdf_georef(ds, data.shape)

        # Should default to EPSG:4326
        assert geo_info["crs"] == "EPSG:4326"


class TestImportErrors:
    """Test error handling for missing dependencies."""

    def test_read_geotiff_no_rasterio(self, monkeypatch):
        """Test read_geotiff error when rasterio not available."""
        # Temporarily disable rasterio
        monkeypatch.setattr(raster, "HAS_RASTERIO", False)

        with pytest.raises(ImportError, match="rasterio required"):
            raster.read_geotiff("dummy.tif")

    def test_write_geotiff_no_rasterio(self, monkeypatch):
        """Test write_geotiff error when rasterio not available."""
        monkeypatch.setattr(raster, "HAS_RASTERIO", False)

        data = np.ones((10, 10))
        with pytest.raises(ImportError, match="rasterio required"):
            raster.write_geotiff(data, "dummy.tif")

    def test_update_netcdf_no_h5netcdf(self, monkeypatch):
        """Test update_netcdf_variable error when h5netcdf not available."""
        monkeypatch.setattr(raster, "HAS_H5NETCDF", False)

        data = np.ones((10, 10))
        with pytest.raises(ImportError, match="h5netcdf required"):
            raster.update_netcdf_variable("dummy.nc", "var", data)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
