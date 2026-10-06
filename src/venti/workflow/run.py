# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Main workflow runner for InSAR displacement calibration.

This module provides workflow execution functions for:
1. Data staging: Downloading and preparing ancillary data for a given OPERA frame
2. Calibration: Calibrating displacement products using GNSS reference data
3. Decomposition: Converting LOS displacement to ENU components
"""

from __future__ import annotations

import logging
from enum import StrEnum
from pathlib import Path
from typing import Literal, cast

from .calibration import CalibrationState, CalibrationWorkflow
from .config import AlgorithmParameters, RunConfig, VentiConfig, WorkflowConfig
from .decomposition import DecompositionState, DecompositionWorkflow

logger = logging.getLogger(__name__)


class WorkflowType(StrEnum):
    """Available workflow types."""

    calibrate = "calibrate"
    decompose = "decompose"


def run_workflow(
    config: WorkflowConfig,
    workflow_type: WorkflowType = WorkflowType.calibrate,
    n_workers: int = 1,
) -> CalibrationState | DecompositionState:
    """Run a workflow using provided configuration.

    Parameters
    ----------
    config : WorkflowConfig
        Configuration object with all workflow parameters
    workflow_type : WorkflowType, optional
        Type of workflow to run: 'calibrate' or 'decompose', by default 'calibrate'
    n_workers : int, optional
        Products calibrated concurrently (see `CalibrationWorkflow.run`), by
        default 1. Not used by the decomposition.

    Returns
    -------
    CalibrationState or DecompositionState
        Final workflow state with processing statistics

    Examples
    --------
    Run calibration workflow::

        from venti.workflow import WorkflowType, load_config, run_workflow

        config = load_config('runconfig.yaml')
        state = run_workflow(config, WorkflowType.calibrate)
        print(f"Processed {state.n_files_processed} files")

    Run decomposition workflow::

        config = load_config('runconfig.yaml')
        state = run_workflow(config, WorkflowType.decompose)

    """
    # Set up worker settings
    import os

    if config.worker_settings.threads_per_worker > 1:
        os.environ["OMP_NUM_THREADS"] = str(config.worker_settings.threads_per_worker)

    # Dispatch to appropriate workflow
    if workflow_type == WorkflowType.calibrate:
        calib_workflow = CalibrationWorkflow(config=config)
        return calib_workflow.run(n_workers=n_workers)
    elif workflow_type == WorkflowType.decompose:
        decomp_workflow = DecompositionWorkflow(config=config)
        return decomp_workflow.run()
    else:
        msg = f"Unknown workflow type: {workflow_type}"
        raise ValueError(msg)


def calibrate_timeseries(
    input_dir: Path,
    los_file: Path,
    mask_file: Path,
    output_dir: Path,
    correction_dir: Path | None = None,
    reference_point: tuple[int, int] | None = None,
    grid_type: Literal["constant", "variable"] = "constant",
    downsample_factor: int = 1,
    window_size_meters: float = 30000,
    posting_meters: float = 30,
    reference_frame: str = "IGS20",
    unwrap_error_correction: bool = True,
    apply_tropo_correction: bool = True,
    apply_solid_earth_tide_correction: bool = True,
    n_workers: int = 1,
    keep_paths_relative: bool = False,
    log_file: Path | None = None,
) -> CalibrationState:
    """Calibrate InSAR displacement timeseries using GNSS reference data.

    Builds a configuration from keyword arguments and runs `CalibrationWorkflow`;
    use `load_config` and `run_workflow` to run from YAML files.

    Parameters
    ----------
    input_dir : Path
        Directory containing input NetCDF displacement files
    los_file : Path
        Path to LOS unit vector file (3-band GeoTIFF: east, north, up)
    mask_file : Path
        Path to mask file (1=valid, 0=invalid)
    output_dir : Path
        Output directory for corrected displacements and intermediate products
    correction_dir : Path, optional
        Directory with tropospheric correction files, by default None
    reference_point : tuple of int, optional
        Reference point (row, col), by default None (auto-select)
    grid_type : str, optional
        GNSS grid type: "constant" (velocity) or "variable" (epoch-specific),
        by default "constant"
    downsample_factor : int, optional
        Downsampling factor for faster processing, by default 1 (no downsampling)
    window_size_meters : float, optional
        Window size for surface fitting in meters, by default 30000
    posting_meters : float, optional
        Pixel posting in meters, by default 30
    reference_frame : str, optional
        GNSS reference frame, by default "IGS20"
    unwrap_error_correction : bool, optional
        Whether to correct islands for unwrap errors, by default True
    apply_tropo_correction : bool, optional
        Whether to apply tropospheric correction when `correction_dir` is
        given, by default True. Set False to skip tropo without unsetting
        `correction_dir`, e.g. for an A/B comparison.
    apply_solid_earth_tide_correction : bool, optional
        Whether to remove the product's ``/corrections/solid_earth_tide``
        layer before the fit, by default True.
    n_workers : int, optional
        Products calibrated concurrently (see `CalibrationWorkflow.run`),
        by default 1.
    keep_paths_relative : bool, optional
        Don't resolve filepaths that are given as relative to be absolute,
        by default False
    log_file : Path, optional
        Also write log messages to this file, by default None.

    Returns
    -------
    CalibrationState
        Final workflow state with processing statistics

    Examples
    --------
    Calibrate with constant velocity model::

        state = calibrate_timeseries(
            input_dir=Path('displacement/'),
            los_file=Path('los_vectors.tif'),
            mask_file=Path('water_mask.tif'),
            output_dir=Path('calibrated/'),
            grid_type='constant'
        )

    Calibrate with tropospheric corrections::

        state = calibrate_timeseries(
            input_dir=Path('displacement/'),
            los_file=Path('los_vectors.tif'),
            mask_file=Path('water_mask.tif'),
            output_dir=Path('calibrated/'),
            correction_dir=Path('tropo_corrections/'),
            grid_type='constant'
        )

    Notes
    -----
    One calibration surface GeoTIFF is written per product, named
    ``<product>_calibration_surface_<grid_type>_<frame>[_downsampleN][_tropo]
    [_set][_nowrap].tif``; subtract it from the raw ``/displacement``.

    """
    from .config import (
        CalibrationInputGroup,
        CalibrationOptions,
        PrimaryExecutable,
        ProcessingOptions,
        ProductPathGroup,
    )

    run_config = RunConfig(
        calibration_input_group=CalibrationInputGroup(
            input_files=input_dir,
            los_file=los_file,
            water_mask=mask_file,
            tropo_files=correction_dir,
            reference_point=reference_point,
        ),
        product_path_group=ProductPathGroup(
            product_path=output_dir,
            scratch_path=output_dir / "scratch",
            sas_output_path=output_dir,
            product_version="1.0",
        ),
        primary_executable=PrimaryExecutable(
            product_type="CAL",
            workflow_name="calibrate",
        ),
        keep_paths_relative=keep_paths_relative,
        log_file=str(log_file) if log_file is not None else None,
    )

    algorithm_params = AlgorithmParameters(
        processing_options=ProcessingOptions(
            cal_downsample_factor=downsample_factor,
        ),
        calibration_options=CalibrationOptions(
            grid_type=grid_type,
            reference_frame=reference_frame,
            unwrap_error_correction=unwrap_error_correction,
            apply_tropo_correction=apply_tropo_correction,
            apply_solid_earth_tide_correction=apply_solid_earth_tide_correction,
            window_size_meters=window_size_meters,
            posting_meters=posting_meters,
        ),
    )

    config = VentiConfig(
        run_config=run_config,
        algorithm_parameters=algorithm_params,
    )

    return cast(
        CalibrationState,
        run_workflow(config, WorkflowType.calibrate, n_workers=n_workers),
    )


def decompose_timeseries(
    input_dir: Path,
    los_file: Path,
    mask_file: Path,
    output_dir: Path,
    asc_dir: Path | None = None,
    desc_dir: Path | None = None,
    reference_point: tuple[int, int] | None = None,
    gpu_enabled: bool = False,
    block_shape: list[int] | None = None,
    keep_paths_relative: bool = False,
    log_file: Path | None = None,
) -> DecompositionState:
    """Decompose InSAR LOS displacement to East-North-Up components.

    Builds a configuration from keyword arguments and runs `DecompositionWorkflow`;
    use `load_config` and `run_workflow` to run from YAML files.

    Parameters
    ----------
    input_dir : Path
        Directory containing input NetCDF displacement files
    los_file : Path
        Path to LOS unit vector file (3-band GeoTIFF: east, north, up)
    mask_file : Path
        Path to mask file (1=valid, 0=invalid)
    output_dir : Path
        Output directory for decomposed displacement components
    asc_dir : Path, optional
        Directory with ascending geometry displacement files, by default None
    desc_dir : Path, optional
        Directory with descending geometry displacement files, by default None
    reference_point : tuple of int, optional
        Reference point (row, col), by default None (auto-select)
    gpu_enabled : bool, optional
        Whether to use GPU for processing (if available), by default False
    block_shape : list of int, optional
        Size (rows, columns) of blocks of data to load at a time, by default [512, 512]
    keep_paths_relative : bool, optional
        Don't resolve filepaths that are given as relative to be absolute,
        by default False
    log_file : Path, optional
        Also write log messages to this file, by default None.

    Returns
    -------
    DecompositionState
        Final workflow state with processing statistics

    Examples
    --------
    Decompose with single geometry::

        state = decompose_timeseries(
            input_dir=Path('displacement/'),
            los_file=Path('los_vectors.tif'),
            mask_file=Path('water_mask.tif'),
            output_dir=Path('decomposed/')
        )

    Decompose with ascending and descending::

        state = decompose_timeseries(
            input_dir=Path('displacement/'),
            los_file=Path('los_vectors.tif'),
            mask_file=Path('water_mask.tif'),
            output_dir=Path('decomposed/'),
            asc_dir=Path('ascending/'),
            desc_dir=Path('descending/')
        )

    Notes
    -----
    This is a PLACEHOLDER implementation. Full decomposition requires:
    - Multiple viewing geometries (ascending/descending)
    - Proper least-squares inversion
    - Uncertainty quantification
    - Quality metrics

    """
    from .config import (
        DecompositionInputGroup,
        PrimaryExecutable,
        ProductPathGroup,
        WorkerSettings,
    )

    # Set default block_shape if not provided
    if block_shape is None:
        block_shape = [512, 512]

    # Build configuration from parameters
    # Note: This API is simplified and assumes asc/desc dirs contain both
    # displacement files and LOS files with standard naming
    run_config = RunConfig(
        decomposition_input_group=DecompositionInputGroup(
            asc_displacement_files=asc_dir if asc_dir else input_dir,
            desc_displacement_files=desc_dir if desc_dir else input_dir,
            asc_los_file=los_file,  # Simplified - should be asc-specific
            desc_los_file=los_file,  # Simplified - should be desc-specific
            water_mask=mask_file,
            reference_point=reference_point,
        ),
        product_path_group=ProductPathGroup(
            product_path=output_dir,
            scratch_path=output_dir / "scratch",
            sas_output_path=output_dir,
            product_version="1.0",
        ),
        primary_executable=PrimaryExecutable(
            product_type="VLM",
            workflow_name="decompose",
        ),
        worker_settings=WorkerSettings(
            gpu_enabled=gpu_enabled,
            threads_per_worker=1,
            block_shape=block_shape,
        ),
        keep_paths_relative=keep_paths_relative,
        log_file=str(log_file) if log_file is not None else None,
    )

    algorithm_params = AlgorithmParameters()

    config = VentiConfig(
        run_config=run_config,
        algorithm_parameters=algorithm_params,
    )

    # Run workflow using the DecompositionWorkflow class
    return cast(
        DecompositionState, run_workflow(config, workflow_type=WorkflowType.decompose)
    )


def run_data_staging(
    frame_id: int,
    date: str,
    output_dir: Path,
    num_workers: int = 4,
    dem_buffer: float = 10_000.0,
    skip_tropo: bool = False,
    skip_gnss: bool = False,
    gnss_reference_frame: str = "IGS20",
    gnss_padding: float = 0.0,
    gnss_start_year: float = 2014.0,
    log_file: Path | None = None,
) -> None:
    """Stage all ancillary data for a single OPERA DISP-S1 frame.

    Downloads the DISP-S1 product and all ancillary data needed for downstream
    calibration: DEM, LOS geometry, tropospheric corrections, and UNR GNSS
    velocities. DEM and LOS outputs are idempotent — existing files are reused.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    date : str
        Secondary date of the interferogram to stage (YYYY-MM-DD or YYYYMMDD).
    output_dir : Path
        Root directory for all staged outputs.
    num_workers : int, optional
        Number of parallel workers for downloads and processing. Default is 4.
    dem_buffer : float, optional
        Buffer in meters around the frame extent for DEM generation.
        Default is 10,000 m (10 km).
    skip_tropo : bool, optional
        Skip tropospheric correction processing. Default is False.
    skip_gnss : bool, optional
        Skip UNR GNSS download and velocity estimation. Default is False.
    gnss_reference_frame : str, optional
        GNSS reference frame for UNR data (``'IGS20'`` or ``'IGS14'``).
        Default is ``'IGS20'``.
    gnss_padding : float, optional
        Extra padding in meters beyond frame bounds when searching for GNSS
        stations. Default is ``0.0``.
    gnss_start_year : float, optional
        Exclude GNSS observations before this decimal year when estimating
        velocities. Default is ``2014.0``.
    log_file : Path, optional
        Also write log messages to this file, by default None.

    Examples
    --------
    Stage a frame with all ancillary data::

        run_data_staging(
            frame_id=8887,
            date="2016-06-15",
            output_dir=Path("./data"),
        )

    Stage without tropospheric corrections or GNSS::

        run_data_staging(
            frame_id=8887,
            date="2016-06-15",
            output_dir=Path("./data"),
            skip_tropo=True,
            skip_gnss=True,
        )

    """
    from ..log_setup import configure_logging
    from .stage_frame_data import stage_frame

    configure_logging(log_file=log_file)
    stage_frame(
        frame_id=frame_id,
        date=date,
        output_dir=output_dir,
        num_workers=num_workers,
        dem_buffer=dem_buffer,
        skip_tropo=skip_tropo,
        skip_gnss=skip_gnss,
        gnss_reference_frame=gnss_reference_frame,
        gnss_padding=gnss_padding,
        gnss_start_year=gnss_start_year,
    )
    logger.info("Data staging complete for frame %d on %s", frame_id, date)


def run_data_staging_window(
    frame_id: int,
    start: str,
    end: str,
    output_dir: Path,
    num_workers: int = 4,
    dem_buffer: float = 10_000.0,
    skip_tropo: bool = False,
    skip_gnss: bool = False,
    gnss_reference_frame: str = "IGS20",
    gnss_padding: float = 0.0,
    gnss_start_year: float = 2014.0,
    log_file: Path | None = None,
) -> list[Path]:
    """Stage all ancillary data for multiple DISP-S1 products in a time window.

    Downloads every DISP-S1 product whose secondary date falls within
    ``[start, end]`` and generates all ancillary data required for downstream
    calibration.  Frame-level assets (DEM, LOS, GNSS) are produced once and
    shared across all products; tropospheric corrections are batched over all
    unique epoch sensing times and combined into per-product differential files.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    start : str
        Start of the secondary-date window (YYYY-MM-DD or YYYYMMDD).
    end : str
        End of the secondary-date window (YYYY-MM-DD or YYYYMMDD).
    output_dir : Path
        Root directory for all staged outputs.
    num_workers : int, optional
        Number of parallel workers for downloads and processing. Default is 4.
    dem_buffer : float, optional
        Buffer in metres around the frame extent for DEM generation.
        Default is 10,000 m (10 km).
    skip_tropo : bool, optional
        Skip tropospheric correction processing. Default is False.
    skip_gnss : bool, optional
        Skip UNR GNSS download and velocity estimation. Default is False.
    gnss_reference_frame : str, optional
        GNSS reference frame for UNR data (``'IGS20'`` or ``'IGS14'``).
        Default is ``'IGS20'``.
    gnss_padding : float, optional
        Extra padding in metres beyond frame bounds when searching for GNSS
        stations. Default is ``0.0``.
    gnss_start_year : float, optional
        Exclude GNSS observations before this decimal year when estimating
        velocities. Default is ``2014.0``.
    log_file : Path, optional
        Also write log messages to this file, by default None.

    Returns
    -------
    list[Path]
        Paths of the downloaded DISP-S1 NetCDF files, sorted by secondary date.

    Examples
    --------
    Stage all products for a frame in a given year::

        run_data_staging_window(
            frame_id=8887,
            start="2016-06-01",
            end="2016-12-31",
            output_dir=Path("./data"),
        )

    """
    from ..log_setup import configure_logging
    from .stage_frame_data import stage_window

    configure_logging(log_file=log_file)
    disp_files = stage_window(
        frame_id=frame_id,
        start=start,
        end=end,
        output_dir=output_dir,
        num_workers=num_workers,
        dem_buffer=dem_buffer,
        skip_tropo=skip_tropo,
        skip_gnss=skip_gnss,
        gnss_reference_frame=gnss_reference_frame,
        gnss_padding=gnss_padding,
        gnss_start_year=gnss_start_year,
    )
    logger.info(
        "Data staging complete: %d products for frame %d (%s to %s)",
        len(disp_files),
        frame_id,
        start,
        end,
    )
    return disp_files


if __name__ == "__main__":
    import tyro

    # Support subcommands: calibrate or decompose
    tyro.extras.subcommand_cli_from_dict(
        {
            WorkflowType.calibrate: calibrate_timeseries,
            WorkflowType.decompose: decompose_timeseries,
        }
    )
