"""Unit tests for unwrap corrections module."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

# Optional rasterio import
try:
    from rasterio.transform import Affine

    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    Affine = None

from venti.unwrap import UnwrapCorrector, correct_region_offset


class TestUnwrapCorrector:
    """Test cases for UnwrapCorrector class."""

    def test_init_default(self):
        """Test initialization with default parameters."""
        corrector = UnwrapCorrector()
        assert corrector.min_region_area == 20
        assert corrector.wavelength == 0.0555

    def test_init_custom(self):
        """Test initialization with custom parameters."""
        corrector = UnwrapCorrector(min_region_area=50, wavelength=0.06)
        assert corrector.min_region_area == 50
        assert corrector.wavelength == 0.06

    def test_prepare_displacement(self):
        """Test displacement preparation and masking."""
        corrector = UnwrapCorrector()

        # Create test data with some zeros and NaNs
        disp = np.array([[1.0, 2.0, 0.0], [3.0, np.nan, 4.0], [5.0, 6.0, 7.0]])

        scaled_disp, disp_mask = corrector._prepare_displacement(disp)

        # Check that the displacement was scaled
        assert scaled_disp is not None
        assert scaled_disp.shape == disp.shape

        # Check that invalid values (NaN) are identified
        # disp_mask can be a scalar False if no invalid values, or an array if there are
        if isinstance(disp_mask, np.ndarray):
            assert disp_mask[1, 1]  # NaN value should be masked

        # Check that the scaled displacement is valid (no NaNs after processing)
        assert not np.isnan(scaled_disp).any()

    def test_watershed_segmentation(self):
        """Test watershed segmentation."""
        corrector = UnwrapCorrector(min_region_area=10)

        # Create synthetic displacement with distinct regions
        disp = np.ones((50, 50)) * 5.0
        disp[10:20, 10:20] = 1.0  # Region 1 (lower values)
        disp[30:40, 30:40] = 10.0  # Region 2 (higher values)

        mask = np.ones((50, 50), dtype=bool)

        # Prepare displacement first
        scaled_disp, _disp_mask = corrector._prepare_displacement(disp)

        labeled, valid_labels = corrector._watershed_segmentation(scaled_disp, mask)

        # Should have at least 1 valid region
        assert len(valid_labels) >= 1
        assert labeled.shape == (50, 50)

    def test_compute_regional_medians(self):
        """Test regional median computation."""
        corrector = UnwrapCorrector()

        # Create labeled regions
        labeled = np.array([[1, 1, 2, 2], [1, 1, 2, 2], [3, 3, 3, 3]])

        disp = np.array(
            [[1.0, 1.0, 5.0, 5.0], [1.0, 1.0, 5.0, 5.0], [10.0, 10.0, 10.0, 10.0]]
        )

        valid_labels = np.array([1, 2, 3])

        medians = corrector._compute_regional_medians(disp, labeled, valid_labels)

        assert len(medians) == 3
        assert medians[0] == 1.0  # Median of region 1
        assert medians[1] == 5.0  # Median of region 2
        assert medians[2] == 10.0  # Median of region 3

    def test_compute_unwrap_cycles(self):
        """Test unwrap cycle computation."""
        corrector = UnwrapCorrector(wavelength=0.0555)

        # Create medians with known offsets
        medians = np.array([0.0, 0.0555, -0.0555, 0.111, 0.02])

        cycles = corrector._compute_unwrap_cycles(medians, ref_label=0)

        assert cycles[0] == 0  # Reference
        assert cycles[1] == 1  # One wavelength above
        assert cycles[2] == -1  # One wavelength below
        assert cycles[3] == 2  # Two wavelengths above
        assert cycles[4] == 0  # Less than half wavelength, rounds to 0

    def test_compute_unwrap_cycles_with_nan(self):
        """Test unwrap cycle computation with NaN values."""
        corrector = UnwrapCorrector(wavelength=0.0555)

        medians = np.array([0.0, 0.0555, np.nan, 0.111])
        cycles = corrector._compute_unwrap_cycles(medians, ref_label=0)

        # NaN should be converted to 0
        assert cycles[2] == 0
        assert not np.isnan(cycles).any()

    def test_correct_simple(self):
        """An island one wavelength off is brought back to the reference island."""
        rng = np.random.default_rng(0)
        disp = np.full((40, 60), np.nan)
        disp[5:35, 3:25] = 0.010 + 0.002 * rng.standard_normal((30, 22))
        # Separated by no-data, and one full wavelength off.
        disp[5:35, 33:57] = 0.010 + 0.0555 + 0.002 * rng.standard_normal((30, 24))

        corrector = UnwrapCorrector(wavelength=0.0555)
        corrected = np.ma.filled(corrector.correct(disp, np.isfinite(disp)), np.nan)

        np.testing.assert_array_equal(corrector.unwrap_cycles_, [0, 1])
        np.testing.assert_allclose(
            corrected[5:35, 33:57], disp[5:35, 33:57] - 0.0555, atol=1e-6
        )
        np.testing.assert_allclose(corrected[5:35, 3:25], disp[5:35, 3:25], atol=1e-6)

    def test_small_islands_keep_their_values(self):
        """Pixels outside corrected regions must not be left uninitialised.

        Regression: islands below `min_region_area` in no-data connected to
        the image edge (e.g. small islands at sea) came out as leftover
        memory instead of their own values.
        """
        rng = np.random.default_rng(0)
        disp = 0.003 + 0.001 * rng.standard_normal((60, 80))
        disp[:, 40:] = np.nan  # open sea reaching the image edge
        islands = [(15, 60), (30, 66), (45, 53)]
        for r, c in islands:
            disp[r : r + 3, c : c + 3] = 0.002  # 9 px < min_region_area

        corrected = np.ma.filled(
            UnwrapCorrector(wavelength=0.0555).correct(disp, np.isfinite(disp)),
            np.nan,
        )

        valid = np.isfinite(disp)
        np.testing.assert_allclose(corrected[valid], disp[valid], atol=1e-6)


class TestCorrectRegionOffset:
    """Test cases for correct_region_offset convenience function."""

    def test_correct_with_arrays(self):
        """Test correction with array inputs."""
        # Create synthetic data
        disp = np.random.rand(50, 50)
        mask = np.ones((50, 50), dtype=bool)
        mask[:5, :5] = False

        corrected = correct_region_offset(
            input_disp=disp, mask=mask, wavelength=0.0555, min_region_area=20
        )

        # A plain array with NaN where masked, usable by other Venti functions.
        assert not isinstance(corrected, np.ma.MaskedArray)
        assert corrected.shape == disp.shape
        assert np.isnan(corrected[:5, :5]).all()
        assert np.isfinite(corrected[mask]).all()

    def test_correct_with_netcdf(self):
        """Test correction with NetCDF file input."""
        # Create temporary NetCDF file
        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
            tmp_path = tmp.name

        try:
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

            ds.to_netcdf(tmp_path)

            # Test correction
            corrected = correct_region_offset(
                input_disp=tmp_path, wavelength=0.0555, min_region_area=20
            )

            assert corrected is not None
            assert corrected.shape == disp_data.shape

        finally:
            Path(tmp_path).unlink()

    @pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
    def test_correct_with_output_file(self):
        """Test correction with GeoTIFF output."""
        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp_nc:
            nc_path = tmp_nc.name
        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_tif:
            tif_path = tmp_tif.name

        try:
            # Create synthetic NetCDF with georeferencing
            x = np.arange(0, 100, 10)
            y = np.arange(0, 50, 10)
            disp_data = np.random.rand(len(y), len(x))
            mask_data = np.ones((len(y), len(x)), dtype=bool)

            ds = xr.Dataset(
                {
                    "displacement": (["y", "x"], disp_data),
                    "water_mask": (["y", "x"], mask_data),
                    "spatial_ref": ([], 0),
                },
                coords={"x": x, "y": y},
            )

            ds["spatial_ref"].attrs["GeoTransform"] = "0.0 10.0 0.0 50.0 0.0 -10.0"
            ds["spatial_ref"].attrs["crs_wkt"] = (
                'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS'
                ' 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
            )

            ds.to_netcdf(nc_path)

            # Test correction with output
            corrected = correct_region_offset(
                input_disp=nc_path,
                wavelength=0.0555,
                min_region_area=20,
                output_file=tif_path,
            )

            assert corrected is not None
            assert Path(tif_path).exists()

        finally:
            Path(nc_path).unlink()
            if Path(tif_path).exists():
                Path(tif_path).unlink()


@pytest.mark.skipif(not HAS_RASTERIO, reason="rasterio not installed")
class TestSaveGeoTIFF:
    """Test cases for save_geotiff method."""

    def test_save_geotiff_with_transform_and_crs(self):
        """Test saving GeoTIFF with explicit transform and CRS."""
        corrector = UnwrapCorrector()

        # Create synthetic corrected data
        corrected = np.ma.masked_array(
            np.random.rand(50, 50), mask=np.zeros((50, 50), dtype=bool)
        )

        transform = Affine(10.0, 0.0, 0.0, 0.0, -10.0, 500.0)
        crs = "EPSG:32611"

        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            corrector.save_geotiff(
                corrected_disp=corrected,
                output_path=tmp_path,
                transform=transform,
                crs=crs,
                nodata=-9999,
            )

            assert Path(tmp_path).exists()

        finally:
            Path(tmp_path).unlink()

    def test_save_geotiff_missing_params(self):
        """Test error handling when transform or CRS is missing."""
        corrector = UnwrapCorrector()

        corrected = np.ma.masked_array(
            np.random.rand(50, 50), mask=np.zeros((50, 50), dtype=bool)
        )

        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            # Should raise ValueError when transform is missing
            with pytest.raises(ValueError, match=r"transform|reference"):
                corrector.save_geotiff(
                    corrected_disp=corrected, output_path=tmp_path, crs="EPSG:32611"
                )

        finally:
            if Path(tmp_path).exists():
                Path(tmp_path).unlink()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
