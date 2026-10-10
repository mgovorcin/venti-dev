# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Unit tests for workflow module.

This module tests configuration handling and utility functions
for the Venti workflow system.
"""

import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from venti.spatial.resample import downsample_array, upsample_array
from venti.workflow.config import (
    AlgorithmParameters,
    CalibrationInputGroup,
    CalibrationOptions,
    DecompositionInputGroup,
    DecompositionOptions,
    OutputOptions,
    PrimaryExecutable,
    ProcessingOptions,
    ProductPathGroup,
    RunConfig,
    VentiConfig,
    WorkerSettings,
    create_config_templates,
    load_config,
)
from venti.workflow.utils import (
    datetime_to_decimal_year,
    extract_dates_from_filename,
    match_correction_to_displacement,
    parse_window_size_meters,
)


class TestProcessingOptions:
    """Test cases for ProcessingOptions model."""

    def test_default_values(self):
        """Test default processing options."""
        opts = ProcessingOptions()
        assert opts.cal_downsample_factor == 1
        assert opts.downsample_method == "mean"
        assert opts.downsample_weighted is False
        assert opts.vlm_output_posting_meters == 120.0

    def test_custom_values(self):
        """Test custom processing options."""
        opts = ProcessingOptions(
            cal_downsample_factor=2,
            downsample_method="median",
            downsample_weighted=True,
            vlm_output_posting_meters=90.0,
        )
        assert opts.cal_downsample_factor == 2
        assert opts.downsample_method == "median"
        assert opts.downsample_weighted is True
        assert opts.vlm_output_posting_meters == 90.0

    def test_invalid_method(self):
        """Test validation of downsample_method."""
        with pytest.raises(ValueError, match=r"downsample_method|Input should be"):
            ProcessingOptions(downsample_method="invalid")

    def test_invalid_downsample_factor(self):
        """Test validation of cal_downsample_factor."""
        with pytest.raises(
            ValueError, match=r"greater than or equal to|Input should be"
        ):
            ProcessingOptions(cal_downsample_factor=0)


class TestCalibrationOptions:
    """Test cases for CalibrationOptions model."""

    def test_default_values(self):
        """Test default calibration options."""
        opts = CalibrationOptions()
        assert opts.grid_type == "constant"
        assert opts.reference_frame == "IGS20"
        assert (
            opts.unwrap_error_correction is False
        )  # off until TS-U1 (PRD R-U1), schema v2
        assert opts.window_size_meters == 30000.0
        assert opts.posting_meters == 30.0
        assert opts.apply_tropo_correction is True
        assert opts.apply_solid_earth_tide_correction is True

    def test_custom_values(self):
        """Test custom calibration options."""
        opts = CalibrationOptions(
            grid_type="variable",
            reference_frame="IGS14",
            unwrap_error_correction=False,
            window_size_meters=50000.0,
            posting_meters=60.0,
        )
        assert opts.grid_type == "variable"
        assert opts.reference_frame == "IGS14"
        assert opts.unwrap_error_correction is False
        assert opts.window_size_meters == 50000.0
        assert opts.posting_meters == 60.0

    def test_invalid_grid_type(self):
        """Test validation of grid_type."""
        with pytest.raises(ValueError, match=r"grid_type|Input should be"):
            CalibrationOptions(grid_type="invalid")

    def test_smoothing_method_default(self):
        opts = CalibrationOptions()
        assert opts.calibration_surface_smoothing_method == "gaussian"

    def test_smoothing_method_valid_values(self):
        for method in ("gaussian", "gaussian_fft", "hanning_fft", "savitzky_golay"):
            opts = CalibrationOptions(calibration_surface_smoothing_method=method)
            assert opts.calibration_surface_smoothing_method == method

    def test_smoothing_method_invalid_value(self):
        with pytest.raises(ValueError, match="Input should be"):
            CalibrationOptions(calibration_surface_smoothing_method="box_filter")

    def test_smoothing_sigma_default_is_none(self):
        opts = CalibrationOptions()
        assert opts.calibration_surface_smoothing_sigma is None

    def test_smoothing_sigma_zero_accepted(self):
        """Zero is the sentinel for disabling smoothing — must be accepted."""
        opts = CalibrationOptions(calibration_surface_smoothing_sigma=0)
        assert opts.calibration_surface_smoothing_sigma == 0

    def test_smoothing_sigma_negative_rejected(self):
        with pytest.raises(
            ValueError, match=r"greater than or equal to|Input should be"
        ):
            CalibrationOptions(calibration_surface_smoothing_sigma=-1.0)

    def test_event_mask_buffer_default_is_zero(self):
        opts = CalibrationOptions()
        assert opts.event_mask_buffer_pixels == 0

    def test_event_mask_buffer_custom(self):
        opts = CalibrationOptions(event_mask_buffer_pixels=15)
        assert opts.event_mask_buffer_pixels == 15

    def test_event_mask_buffer_negative_rejected(self):
        with pytest.raises(
            ValueError, match=r"greater than or equal to|Input should be"
        ):
            CalibrationOptions(event_mask_buffer_pixels=-1)


class TestDecompositionOptions:
    """Test cases for DecompositionOptions model."""

    def test_default_values(self):
        """Test default decomposition options."""
        opts = DecompositionOptions()
        assert opts.inversion_method == "weighted_least_squares"
        assert opts.uncertainty_propagation is True
        assert opts.min_geometries == 2
        assert opts.vertical_only is False
        assert opts.geometry_weighting == "temporal_coherence"
        assert opts.quality_threshold == 0.5

    def test_custom_values(self):
        """Test custom decomposition options."""
        opts = DecompositionOptions(
            inversion_method="ridge",
            uncertainty_propagation=False,
            min_geometries=3,
            vertical_only=True,
            geometry_weighting="uniform",
            quality_threshold=0.7,
        )
        assert opts.inversion_method == "ridge"
        assert opts.uncertainty_propagation is False
        assert opts.min_geometries == 3
        assert opts.vertical_only is True
        assert opts.geometry_weighting == "uniform"
        assert opts.quality_threshold == 0.7

    def test_invalid_quality_threshold(self):
        """Test validation of quality_threshold."""
        with pytest.raises(ValueError, match=r"less than or equal to|Input should be"):
            DecompositionOptions(quality_threshold=1.5)
        with pytest.raises(
            ValueError, match=r"greater than or equal to|Input should be"
        ):
            DecompositionOptions(quality_threshold=-0.1)


class TestOutputOptions:
    """Test cases for OutputOptions model."""

    def test_default_values(self):
        """Test default output options."""
        opts = OutputOptions()
        assert len(opts.gtiff_creation_options) > 0
        assert "COMPRESS=lzw" in opts.gtiff_creation_options
        assert opts.add_overviews is False
        assert opts.compression == "lzw"


class TestProductPathGroup:
    """Test cases for ProductPathGroup model."""

    def test_default_values(self):
        """Test default product path group."""
        paths = ProductPathGroup()
        assert paths.product_path == Path("output")
        assert paths.scratch_path == Path("scratch")
        assert paths.sas_output_path == Path("output")
        assert paths.product_version == "1.0"

    def test_custom_paths(self):
        """Test custom product paths."""
        paths = ProductPathGroup(
            product_path=Path("/custom/output"),
            scratch_path=Path("/custom/scratch"),
            sas_output_path=Path("/custom/sas"),
            product_version="2.5",
        )
        assert paths.product_path == Path("/custom/output")
        assert paths.scratch_path == Path("/custom/scratch")
        assert paths.sas_output_path == Path("/custom/sas")
        assert paths.product_version == "2.5"

    def test_string_path_conversion(self):
        """Test conversion of string paths to Path objects."""
        paths = ProductPathGroup(
            product_path="custom/output",
            scratch_path="custom/scratch",
        )
        assert isinstance(paths.product_path, Path)
        assert isinstance(paths.scratch_path, Path)


class TestWorkerSettings:
    """Test cases for WorkerSettings model."""

    def test_default_values(self):
        """Test default worker settings."""
        settings = WorkerSettings()
        assert settings.gpu_enabled is False
        assert settings.threads_per_worker == 1
        assert settings.block_shape == [512, 512]

    def test_custom_values(self):
        """Test custom worker settings."""
        settings = WorkerSettings(
            gpu_enabled=True,
            threads_per_worker=4,
            block_shape=[256, 256],
        )
        assert settings.gpu_enabled is True
        assert settings.threads_per_worker == 4
        assert settings.block_shape == [256, 256]

    def test_invalid_block_shape(self):
        """Test validation of block_shape."""
        with pytest.raises(
            ValueError,
            match=r"exactly 2 elements|min_length|List should have at least",
        ):
            WorkerSettings(block_shape=[512])
        with pytest.raises(
            ValueError,
            match=r"exactly 2 elements|max_length|List should have at most",
        ):
            WorkerSettings(block_shape=[512, 512, 512])
        with pytest.raises(
            ValueError, match=r"dimensions must be positive|greater than"
        ):
            WorkerSettings(block_shape=[0, 512])


class TestPrimaryExecutable:
    """Test cases for PrimaryExecutable model."""

    def test_default_values(self):
        """Test default primary executable."""
        exe = PrimaryExecutable()
        assert exe.product_type == "CAL"
        assert exe.workflow_name == "calibrate"

    def test_decomposition_workflow(self):
        """Test decomposition workflow configuration."""
        exe = PrimaryExecutable(
            product_type="VLM",
            workflow_name="decompose",
        )
        assert exe.product_type == "VLM"
        assert exe.workflow_name == "decompose"


class TestCalibrationInputGroup:
    """Test cases for CalibrationInputGroup model."""

    def test_with_placeholder_paths(self):
        """Test with placeholder paths (for templates)."""
        inputs = CalibrationInputGroup(
            input_files=Path("path/to/displacement/files"),
            los_file=Path("path/to/los.tif"),
            water_mask=Path("path/to/mask.tif"),
        )
        assert inputs.input_files == Path("path/to/displacement/files")
        assert inputs.los_file == Path("path/to/los.tif")
        assert inputs.water_mask == Path("path/to/mask.tif")
        assert inputs.tropo_files is None
        assert inputs.reference_point is None

    def test_with_valid_paths(self):
        """Test with valid file paths."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)

            # Create test files and directories
            input_dir = tmpdir / "displacement"
            input_dir.mkdir()
            los_file = tmpdir / "los.tif"
            los_file.touch()
            mask_file = tmpdir / "mask.tif"
            mask_file.touch()

            inputs = CalibrationInputGroup(
                input_files=input_dir,
                los_file=los_file,
                water_mask=mask_file,
            )
            assert inputs.input_files.exists()
            assert inputs.los_file.exists()
            assert inputs.water_mask.exists()

    def test_with_reference_point(self):
        """Test with reference point specified."""
        inputs = CalibrationInputGroup(
            input_files=Path("path/to/displacement/files"),
            los_file=Path("path/to/los.tif"),
            water_mask=Path("path/to/mask.tif"),
            reference_point=(100, 200),
        )
        assert inputs.reference_point == (100, 200)

    def test_missing_file_validation(self):
        """Test validation of missing files."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            nonexistent = tmpdir / "nonexistent.tif"

            with pytest.raises(ValueError, match="does not exist"):
                CalibrationInputGroup(
                    input_files=nonexistent,
                    los_file=Path("path/to/los.tif"),
                    water_mask=Path("path/to/mask.tif"),
                )


class TestDecompositionInputGroup:
    """Test cases for DecompositionInputGroup model."""

    def test_with_placeholder_paths(self):
        """Test with placeholder paths (for templates)."""
        inputs = DecompositionInputGroup(
            asc_displacement_files=Path("path/to/asc"),
            desc_displacement_files=Path("path/to/desc"),
            asc_los_file=Path("path/to/asc_los.tif"),
            desc_los_file=Path("path/to/desc_los.tif"),
            water_mask=Path("path/to/mask.tif"),
        )
        assert inputs.asc_displacement_files == Path("path/to/asc")
        assert inputs.desc_displacement_files == Path("path/to/desc")
        assert inputs.asc_los_file == Path("path/to/asc_los.tif")
        assert inputs.desc_los_file == Path("path/to/desc_los.tif")
        assert inputs.water_mask == Path("path/to/mask.tif")


class TestRunConfig:
    """Test cases for RunConfig model."""

    def test_calibration_workflow(self):
        """Test run config for calibration workflow."""
        config = RunConfig(
            calibration_input_group=CalibrationInputGroup(
                input_files=Path("path/to/displacement/files"),
                los_file=Path("path/to/los.tif"),
                water_mask=Path("path/to/mask.tif"),
            ),
            primary_executable=PrimaryExecutable(
                product_type="CAL",
                workflow_name="calibrate",
            ),
        )
        assert config.primary_executable.workflow_name == "calibrate"
        assert config.calibration_input_group is not None
        assert config.decomposition_input_group is None

    def test_decomposition_workflow(self):
        """Test run config for decomposition workflow."""
        config = RunConfig(
            decomposition_input_group=DecompositionInputGroup(
                asc_displacement_files=Path("path/to/asc"),
                desc_displacement_files=Path("path/to/desc"),
                asc_los_file=Path("path/to/asc_los.tif"),
                desc_los_file=Path("path/to/desc_los.tif"),
                water_mask=Path("path/to/mask.tif"),
            ),
            primary_executable=PrimaryExecutable(
                product_type="VLM",
                workflow_name="decompose",
            ),
        )
        assert config.primary_executable.workflow_name == "decompose"
        assert config.decomposition_input_group is not None
        assert config.calibration_input_group is None

    def test_missing_input_group_validation(self):
        """Test validation when input group is missing."""
        with pytest.raises(ValueError, match="calibration_input_group is required"):
            RunConfig(
                primary_executable=PrimaryExecutable(
                    workflow_name="calibrate",
                ),
            )

    def test_input_file_group_property(self):
        """Test input_file_group property."""
        config = RunConfig(
            calibration_input_group=CalibrationInputGroup(
                input_files=Path("path/to/displacement/files"),
                los_file=Path("path/to/los.tif"),
                water_mask=Path("path/to/mask.tif"),
            ),
            primary_executable=PrimaryExecutable(workflow_name="calibrate"),
        )
        assert config.input_file_group is config.calibration_input_group

    def test_to_yaml(self):
        """Test saving run config to YAML."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = Path(tmpdir) / "runconfig.yaml"

            config = RunConfig(
                calibration_input_group=CalibrationInputGroup(
                    input_files=Path("path/to/displacement/files"),
                    los_file=Path("path/to/los.tif"),
                    water_mask=Path("path/to/mask.tif"),
                ),
                primary_executable=PrimaryExecutable(workflow_name="calibrate"),
            )
            config.to_yaml(yaml_path)

            assert yaml_path.exists()

            # Verify content
            with open(yaml_path) as f:
                data = yaml.safe_load(f)
            assert data is not None
            assert "calibration_input_group" in data


class TestAlgorithmParameters:
    """Test cases for AlgorithmParameters model."""

    def test_default_values(self):
        """Test default algorithm parameters."""
        params = AlgorithmParameters()
        assert isinstance(params.processing_options, ProcessingOptions)
        assert isinstance(params.calibration_options, CalibrationOptions)
        assert isinstance(params.decomposition_options, DecompositionOptions)
        assert isinstance(params.output_options, OutputOptions)

    def test_to_yaml_and_from_yaml(self):
        """Test saving and loading algorithm parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = Path(tmpdir) / "algorithm_parameters.yaml"

            params = AlgorithmParameters(
                calibration_options=CalibrationOptions(
                    grid_type="variable",
                    reference_frame="IGS14",
                ),
            )
            params.to_yaml(yaml_path)

            assert yaml_path.exists()

            # Load back
            loaded = AlgorithmParameters.from_yaml(yaml_path)
            assert loaded.calibration_options.grid_type == "variable"
            assert loaded.calibration_options.reference_frame == "IGS14"


class TestVentiConfig:
    """Test cases for VentiConfig model."""

    def test_from_yaml_files(self):
        """Test loading VentiConfig from separate YAML files."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            runconfig_path = tmpdir / "runconfig.yaml"
            algorithm_params_path = tmpdir / "algorithm_parameters.yaml"

            # Create run config
            run_config = RunConfig(
                calibration_input_group=CalibrationInputGroup(
                    input_files=Path("path/to/displacement/files"),
                    los_file=Path("path/to/los.tif"),
                    water_mask=Path("path/to/mask.tif"),
                ),
                primary_executable=PrimaryExecutable(workflow_name="calibrate"),
            )
            run_config.to_yaml(runconfig_path)

            # Create algorithm parameters
            algorithm_params = AlgorithmParameters()
            algorithm_params.to_yaml(algorithm_params_path)

            # Load combined config
            config = VentiConfig.from_yaml_files(runconfig_path, algorithm_params_path)

            assert isinstance(config.run_config, RunConfig)
            assert isinstance(config.algorithm_parameters, AlgorithmParameters)

    def test_compatibility_properties(self):
        """Test backward compatibility properties (grid_settings, input_options, etc)."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            runconfig_path = tmpdir / "runconfig.yaml"
            algorithm_params_path = tmpdir / "algorithm_parameters.yaml"

            # Create config with specific values
            run_config = RunConfig(
                calibration_input_group=CalibrationInputGroup(
                    input_files=Path("path/to/displacement/files"),
                    los_file=Path("path/to/los.tif"),
                    water_mask=Path("path/to/mask.tif"),
                ),
                primary_executable=PrimaryExecutable(workflow_name="calibrate"),
            )
            run_config.to_yaml(runconfig_path)

            algorithm_params = AlgorithmParameters(
                processing_options=ProcessingOptions(
                    cal_downsample_factor=2,
                    downsample_method="median",
                    downsample_weighted=True,
                    vlm_output_posting_meters=90.0,
                ),
                calibration_options=CalibrationOptions(
                    grid_type="variable",
                ),
            )
            algorithm_params.to_yaml(algorithm_params_path)

            # Load combined config
            config = VentiConfig.from_yaml_files(runconfig_path, algorithm_params_path)

            # Test compatibility properties
            assert config.grid_settings.grid_type == "variable"
            assert config.grid_settings.downsample_factor == 2
            assert config.grid_settings.downsample_method == "median"
            assert config.grid_settings.downsample_weighted is True
            assert config.grid_settings.output_posting_meters == 90.0

            assert config.input_options is not None
            assert config.worker_settings is not None


class TestConfigHelperFunctions:
    """Test cases for configuration helper functions."""

    def test_create_config_templates(self):
        """Test creating configuration template files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            runconfig_path, algorithm_params_path = create_config_templates(tmpdir)

            assert runconfig_path.exists()
            assert algorithm_params_path.exists()

            # Verify runconfig content
            with open(runconfig_path) as f:
                runconfig_data = yaml.safe_load(f)
            assert "calibration_input_group" in runconfig_data
            assert "decomposition_input_group" in runconfig_data
            assert "primary_executable" in runconfig_data

            # Verify algorithm parameters content
            with open(algorithm_params_path) as f:
                algorithm_data = yaml.safe_load(f)
            assert "processing_options" in algorithm_data
            assert "calibration_options" in algorithm_data
            assert "decomposition_options" in algorithm_data
            assert "output_options" in algorithm_data

    def test_load_config(self):
        """Test load_config helper function."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            runconfig_path = tmpdir / "runconfig.yaml"
            algorithm_params_path = tmpdir / "algorithm_parameters.yaml"

            # Create a minimal valid config (only calibration input group)
            run_config = RunConfig(
                calibration_input_group=CalibrationInputGroup(
                    input_files=Path("path/to/displacement/files"),
                    los_file=Path("path/to/los.tif"),
                    water_mask=Path("path/to/mask.tif"),
                ),
                primary_executable=PrimaryExecutable(workflow_name="calibrate"),
            )
            run_config.to_yaml(runconfig_path)

            # Create algorithm parameters
            algorithm_params = AlgorithmParameters()
            algorithm_params.to_yaml(algorithm_params_path)

            # Load config
            config = load_config(runconfig_path, algorithm_params_path)

            assert isinstance(config, VentiConfig)
            assert isinstance(config.run_config, RunConfig)
            assert isinstance(config.algorithm_parameters, AlgorithmParameters)

    def test_load_config_auto_find_algorithm_params(self):
        """Test load_config finding algorithm_parameters.yaml automatically."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            runconfig_path = tmpdir / "runconfig.yaml"
            algorithm_params_path = tmpdir / "algorithm_parameters.yaml"

            # Create a minimal valid config (only calibration input group)
            run_config = RunConfig(
                calibration_input_group=CalibrationInputGroup(
                    input_files=Path("path/to/displacement/files"),
                    los_file=Path("path/to/los.tif"),
                    water_mask=Path("path/to/mask.tif"),
                ),
                primary_executable=PrimaryExecutable(workflow_name="calibrate"),
            )
            run_config.to_yaml(runconfig_path)

            # Create algorithm parameters in same directory
            algorithm_params = AlgorithmParameters()
            algorithm_params.to_yaml(algorithm_params_path)

            # Load config without specifying algorithm_params_path
            config = load_config(runconfig_path)

            assert isinstance(config, VentiConfig)


# ============================================================================
# Workflow Utils Tests
# ============================================================================


class TestExtractDatesFromFilename:
    """Test cases for extract_dates_from_filename function."""

    def test_opera_filename_two_dates(self):
        """Test extracting two dates from OPERA filename."""
        filename = "OPERA_L3_DISP-S1_20200101T000000_20200115T000000.nc"
        ref_date, sec_date = extract_dates_from_filename(filename)

        assert ref_date == datetime(2020, 1, 1).date()
        assert sec_date == datetime(2020, 1, 15).date()

    def test_filename_one_date(self):
        """Test extracting one date from filename."""
        filename = "data_20200115T120000.nc"
        ref_date, sec_date = extract_dates_from_filename(filename)

        assert ref_date is None
        assert sec_date == datetime(2020, 1, 15, 12, 0, 0).date()

    def test_filename_no_dates(self):
        """Test filename with no dates."""
        filename = "data.nc"
        ref_date, sec_date = extract_dates_from_filename(filename)

        assert ref_date is None
        assert sec_date is None

    def test_path_object_input(self):
        """Test with Path object input."""
        filename = Path("OPERA_L3_DISP-S1_20200101T000000_20200115T000000.nc")
        ref_date, sec_date = extract_dates_from_filename(filename)

        assert ref_date == datetime(2020, 1, 1).date()
        assert sec_date == datetime(2020, 1, 15).date()


class TestDatetimeToDecimalYear:
    """Test cases for datetime_to_decimal_year function."""

    def test_start_of_year(self):
        """Test conversion at start of year."""
        dt = datetime(2020, 1, 1, 0, 0, 0)
        decimal_year = datetime_to_decimal_year(dt)
        assert decimal_year == 2020.0

    def test_middle_of_year(self):
        """Test conversion at approximately middle of year."""
        dt = datetime(2020, 7, 1, 0, 0, 0)
        decimal_year = datetime_to_decimal_year(dt)
        # Should be approximately 2020.5 (within rounding)
        assert 2020.49 < decimal_year < 2020.51

    def test_end_of_year(self):
        """Test conversion at end of year."""
        dt = datetime(2020, 12, 31, 23, 59, 59)
        decimal_year = datetime_to_decimal_year(dt)
        # Should be very close to 2021
        assert decimal_year > 2020.99


class TestMatchCorrectionToDisplacement:
    """Test cases for match_correction_to_displacement function."""

    def test_no_corrections(self):
        """Test matching with no correction files."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)

            # Create displacement files
            disp1 = tmpdir / "disp_20200101T000000_20200115T000000.nc"
            disp2 = tmpdir / "disp_20200101T000000_20200130T000000.nc"
            disp1.touch()
            disp2.touch()

            disp_files = [disp1, disp2]
            matches = match_correction_to_displacement(None, disp_files)

            assert len(matches) == 2
            assert all(ref is None and sec is None for ref, sec, _ in matches)

    def test_with_matching_corrections(self):
        """Test matching per-epoch correction files to displacement files."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)

            # Create displacement files (shared reference date, two secondary dates)
            disp1 = tmpdir / "disp_20200101T000000_20200115T000000.nc"
            disp2 = tmpdir / "disp_20200101T000000_20200130T000000.nc"
            disp1.touch()
            disp2.touch()

            # Per-epoch correction files: one per sensing date
            corr_ref = tmpdir / "corr_20200101T000000.tif"  # shared reference epoch
            corr_sec1 = tmpdir / "corr_20200115T000000.tif"  # secondary of disp1
            corr_sec2 = tmpdir / "corr_20200130T000000.tif"  # secondary of disp2
            corr_ref.touch()
            corr_sec1.touch()
            corr_sec2.touch()

            disp_files = [disp1, disp2]
            corr_files = [corr_ref, corr_sec1, corr_sec2]
            matches = match_correction_to_displacement(corr_files, disp_files)

            assert len(matches) == 2
            matches_by_disp = {disp: (ref, sec) for ref, sec, disp in matches}
            assert matches_by_disp[disp1] == (corr_ref, corr_sec1)
            assert matches_by_disp[disp2] == (corr_ref, corr_sec2)


class TestDownsampleArray:
    """Test cases for downsample_array function."""

    def test_downsample_by_factor_2(self):
        """Test downsampling by factor of 2."""
        array = np.ones((100, 100))
        downsampled = downsample_array(array, factor=2)

        assert downsampled.shape == (50, 50)

    def test_downsample_by_factor_1(self):
        """Test downsampling by factor of 1 (no change)."""
        array = np.ones((100, 100))
        downsampled = downsample_array(array, factor=1)

        assert downsampled.shape == array.shape
        assert np.array_equal(downsampled, array)

    def test_downsample_with_nan(self):
        """Test downsampling preserves NaN values."""
        array = np.ones((100, 100))
        array[10:20, 10:20] = np.nan
        downsampled = downsample_array(array, factor=2)

        # Check that NaNs are present in downsampled array
        assert np.any(np.isnan(downsampled))

    def test_downsample_method_mean(self):
        """Test downsampling with mean method."""
        array = np.arange(100).reshape(10, 10).astype(float)
        downsampled = downsample_array(array, factor=2, method="mean")

        assert downsampled.shape == (5, 5)
        # Check that values are reasonable (mean of blocks)
        assert np.all(downsampled >= 0)

    def test_downsample_method_median(self):
        """Test downsampling with median method."""
        array = np.arange(100).reshape(10, 10).astype(float)
        downsampled = downsample_array(array, factor=2, method="median")

        assert downsampled.shape == (5, 5)
        # Check that values are reasonable (median of blocks)
        assert np.all(downsampled >= 0)

    def test_downsample_with_weights(self):
        """Test weighted downsampling."""
        array = np.ones((100, 100))
        weights = np.ones((100, 100)) * 0.5
        # Higher weight in one region
        weights[0:50, 0:50] = 1.0

        downsampled = downsample_array(array, factor=2, method="mean", weights=weights)

        assert downsampled.shape == (50, 50)
        # Should still be close to 1 since data is uniform
        assert np.allclose(downsampled[~np.isnan(downsampled)], 1.0)

    def test_downsample_invalid_method(self):
        """Test that invalid method raises error."""
        array = np.ones((100, 100))
        with pytest.raises(ValueError, match="Invalid method"):
            downsample_array(array, factor=2, method="invalid")


class TestUpsampleArray:
    """Test cases for upsample_array function."""

    def test_upsample_to_target_shape(self):
        """Test upsampling to target shape."""
        array = np.ones((50, 50))
        upsampled = upsample_array(array, target_shape=(100, 100))

        assert upsampled.shape == (100, 100)

    def test_upsample_same_shape(self):
        """Test upsampling to same shape (no change)."""
        array = np.ones((100, 100))
        upsampled = upsample_array(array, target_shape=(100, 100))

        assert upsampled.shape == array.shape
        assert np.array_equal(upsampled, array)

    def test_upsample_with_nan(self):
        """Test upsampling preserves NaN values."""
        array = np.ones((50, 50))
        array[10:20, 10:20] = np.nan
        upsampled = upsample_array(array, target_shape=(100, 100))

        # Check that NaNs are present in upsampled array
        assert np.any(np.isnan(upsampled))


class TestParseWindowSizeMeters:
    """Test cases for parse_window_size_meters function."""

    def test_default_posting(self):
        """Test conversion with default posting (30m)."""
        window_pixels = parse_window_size_meters(30000)
        assert window_pixels == 1000

    def test_custom_posting(self):
        """Test conversion with custom posting."""
        window_pixels = parse_window_size_meters(30000, posting_meters=60)
        assert window_pixels == 500

    def test_rounding(self):
        """Test rounding to nearest pixel."""
        window_pixels = parse_window_size_meters(30001, posting_meters=30)
        assert window_pixels == 1000
