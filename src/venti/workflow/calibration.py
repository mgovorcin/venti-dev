"""Calibration workflow for InSAR displacement products.

This module provides standalone calibration step functions (GNSS setup,
LOS/mask loading, reference point selection, event mask lookup, and
per-file displacement calibration) that take their inputs as explicit
arguments, plus a `CalibrationWorkflow` class that orchestrates them to
run a single product (`run_single`) or a full batch (`run`).
"""

from __future__ import annotations

import logging
import os
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from tqdm import tqdm

from .config import WorkflowConfig
from .utils import with_scratch_temp_dir

if TYPE_CHECKING:
    from ..gnss.reference import GNSSReference
    from ..io.read import RasterReader
    from ..io.write import RasterWriter
    from ..spatial.processor import SpatialProcessor

logger = logging.getLogger(__name__)


@dataclass
class CalibrationState:
    """State tracking for calibration workflow.

    Attributes
    ----------
    n_files_total : int
        Total number of files to process
    n_files_processed : int
        Number of files processed
    n_files_failed : int
        Number of files that failed
    output_files : list[Path]
        List of output file paths

    """

    n_files_total: int = 0
    n_files_processed: int = 0
    n_files_failed: int = 0
    output_files: list[Path] = field(default_factory=list)

    @property
    def success_rate(self) -> float:
        """Calculate success rate."""
        if self.n_files_processed == 0:
            return 0.0
        return (
            100.0
            * (self.n_files_processed - self.n_files_failed)
            / self.n_files_processed
        )


def setup_gnss_reference(
    io_reader: RasterReader,
    input_files: Path,
    product_path: Path,
    reference_frame: str = "IGS20",
    grid_type: str = "constant",
) -> GNSSReference:
    """Set up a GNSS reference dataset for a directory of displacement files.

    Parameters
    ----------
    io_reader : RasterReader
        Reader used to determine the spatial bounds and CRS of the first
        displacement file.
    input_files : Path
        Directory containing the input NetCDF displacement files.
    product_path : Path
        Output product directory; downloaded GNSS station files are written
        to ``product_path / "GNSS"``.
    reference_frame : str, optional
        GNSS reference frame, ``'IGS20'`` or ``'IGS14'``, by default ``'IGS20'``.
    grid_type : str, optional
        UNR grid product to download: ``'constant'`` (precomputed linear
        rates, IGS20 only) or ``'variable'`` (per-epoch positions), by
        default ``'constant'``. Determines which downstream computation in
        `~venti.gnss.reference.compute_gnss_los` is valid for the returned
        `GNSSReference`.

    Returns
    -------
    GNSSReference
        Initialised reference object with station data downloaded.

    Raises
    ------
    FileNotFoundError
        If `input_files` contains no ``.nc`` files.
    ValueError
        If no GNSS stations are found within the displacement file bounds,
        or if `grid_type` is not available for `reference_frame`.

    """
    from ..gnss.reference import GNSSReference

    disp_files = sorted(input_files.glob("*.nc"))
    if not disp_files:
        msg = f"No NetCDF files in {input_files}"
        raise FileNotFoundError(msg)

    bounds = io_reader.get_bounds(disp_files[0], as_latlon=False)

    from pyproj import CRS

    netcdf_data = io_reader.read_netcdf(disp_files[0])
    utm_epsg = CRS.from_user_input(netcdf_data.crs).to_epsg()

    gnss_dir = product_path / "GNSS"
    gnss_ref = GNSSReference(
        bounds=bounds,
        output_dir=gnss_dir,
        reference_frame=reference_frame,
        utm_epsg=utm_epsg,
        grid_type=grid_type,
    )

    n_stations = gnss_ref.download_stations()
    if n_stations == 0:
        msg = "No GNSS stations found in area"
        raise ValueError(msg)

    logger.info(f"GNSS setup complete: {n_stations} stations")
    return gnss_ref


def load_los_and_mask(
    io_reader: RasterReader,
    los_file: Path,
    water_mask: Path,
    custom_mask: Path | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load LOS unit-vector components and a valid-pixel mask.

    Parameters
    ----------
    io_reader : RasterReader
        Reader used to load the GeoTIFF inputs.
    los_file : Path
        3-band GeoTIFF containing the east, north, and up LOS unit-vector
        components.
    water_mask : Path
        GeoTIFF water mask (nonzero = valid land pixel).
    custom_mask : Path, optional
        Additional GeoTIFF mask combined with `water_mask` via logical AND.

    Returns
    -------
    los_east : np.ndarray
        LOS east component.
    los_north : np.ndarray
        LOS north component.
    los_up : np.ndarray
        LOS up component.
    mask : np.ndarray
        Boolean valid-pixel mask.

    Raises
    ------
    ValueError
        If `los_file` is not a 3-band raster.

    """
    logger.info("Loading LOS unit vectors and mask...")

    los_data = io_reader.read_geotiff(los_file)
    if los_data.data.ndim == 3:
        los_east = los_data.data[0]
        los_north = los_data.data[1]
        los_up = los_data.data[2]
    else:
        msg = "LOS file must be 3-band (east, north, up)"
        raise ValueError(msg)

    mask_data = io_reader.read_geotiff(water_mask)
    mask = mask_data.data.astype(bool)

    # Combine with custom mask if provided (logical AND — a pixel must be
    # valid in both masks to be included in calibration)
    if custom_mask is not None:
        custom_data = io_reader.read_geotiff(custom_mask)
        custom = custom_data.data.astype(bool)
        assert custom.shape == mask.shape, (
            f"custom_mask shape {custom.shape} does not match "
            f"water_mask shape {mask.shape}"
        )
        mask = mask & custom
        n_removed = int((~custom & mask_data.data.astype(bool)).sum())
        logger.info("Custom mask applied: %d additional pixels masked", n_removed)

    logger.info("LOS and mask loaded")
    return los_east, los_north, los_up, mask


def find_reference_point(
    input_files: Path,
    product_path: Path,
    mask: np.ndarray,
    reference_point: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Find or use a configured calibration reference point.

    Parameters
    ----------
    input_files : Path
        Directory containing the input NetCDF displacement files, used to
        compute the average temporal coherence for auto-selection.
    product_path : Path
        Output product directory; the intermediate average coherence raster
        is written here.
    mask : np.ndarray
        Boolean valid-pixel mask (True = valid). Masked pixels are excluded
        from auto-selection by zeroing their coherence before scoring.
    reference_point : tuple of int, optional
        If given, returned unchanged. Otherwise a reference point is
        auto-selected from the average temporal coherence.

    Returns
    -------
    tuple of int
        ``(row, col)`` reference point.

    """
    if reference_point is not None:
        logger.info(f"Using configured reference point: {reference_point}")
        return reference_point

    import tempfile

    import rasterio
    from opera_utils.disp import rebase_reference

    from .utils import compute_average_temporal_coherence

    disp_files = sorted(input_files.glob("*.nc"))
    coherence_file = compute_average_temporal_coherence(
        disp_files,
        product_path,
        variable="temporal_coherence",
    )

    # Zero out invalid pixels so they cannot be selected as reference
    with rasterio.open(coherence_file) as src:
        profile = src.profile.copy()
        coherence = src.read(1)

    valid_mask = mask.squeeze().astype(bool)
    assert (
        coherence.shape == valid_mask.shape
    ), f"Coherence shape {coherence.shape} does not match mask shape {valid_mask.shape}"
    coherence_masked = np.where(valid_mask, coherence, 0.0).astype(profile["dtype"])

    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
        masked_path = Path(tmp.name)

    try:
        with rasterio.open(masked_path, "w", **profile) as dst:
            dst.write(coherence_masked, 1)
        ref_point = rebase_reference.find_reference_point(masked_path)
    finally:
        masked_path.unlink(missing_ok=True)

    logger.info(f"Auto-selected reference point: {ref_point}")
    return ref_point


def find_event_mask_file(event_mask_dir: Path | None, disp_file: Path) -> Path | None:
    """Find the per-epoch event mask file matching a displacement file.

    Looks for a GeoTIFF in `event_mask_dir` whose filename starts with the
    displacement file stem, following the naming convention produced by
    ``generate_event_mask.py``: ``<disp_stem>_<geojson_stem>_mask.tif``.

    Parameters
    ----------
    event_mask_dir : Path or None
        Directory containing per-epoch event mask GeoTIFFs, or `None` if no
        event masking is configured.
    disp_file : Path
        Displacement NetCDF file for which to find a matching event mask.

    Returns
    -------
    Path or None
        Path to the matching event mask GeoTIFF, or `None` if `event_mask_dir`
        is `None` or no match is found.

    """
    if event_mask_dir is None:
        return None

    matches = sorted(event_mask_dir.glob(f"{disp_file.stem}_*mask.tif"))
    if not matches:
        logger.debug(f"No event mask found for {disp_file.name}")
        return None

    if len(matches) > 1:
        logger.warning(
            f"Multiple event masks found for {disp_file.name}; using {matches[0].name}"
        )
    return matches[0]


def process_displacement_file(
    config: WorkflowConfig,
    io_reader: RasterReader,
    io_writer: RasterWriter,
    spatial_processor: SpatialProcessor,
    gnss_manager: GNSSReference,
    disp_file: Path,
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    mask: np.ndarray,
    ref_point: tuple[int, int],
    window_size_x: int,
    window_size_y: int,
    wavelength_m: float,
    tropo_ref_file: Path | None = None,
    tropo_sec_file: Path | None = None,
    event_mask_file: Path | None = None,
    fit_n_jobs: int = -1,
    fit_lock: threading.Lock | None = None,
    recompute_gnss: bool = False,
) -> Path | None:
    """Calibrate a single displacement file against a GNSS LOS reference.

    Reads the product, its tropospheric and solid Earth tide corrections and
    the GNSS LOS field, runs `venti.surface.estimate_calibration_surface`
    (the array-only core), and writes the surface as a GeoTIFF. The surface
    is what to subtract from the raw ``/displacement``.

    Parameters
    ----------
    config : WorkflowConfig
        Workflow configuration.
    io_reader : RasterReader
        Reader for the displacement, tropospheric, solid Earth tide, and
        event mask rasters.
    io_writer : RasterWriter
        Writer used to save the calibration surface GeoTIFF.
    spatial_processor : SpatialProcessor
        Processor used to fit the windowed calibration surface.
    gnss_manager : GNSSReference
        Initialised GNSS reference data with stations already downloaded.
    disp_file : Path
        Displacement file path.
    los_east : np.ndarray
        LOS east component.
    los_north : np.ndarray
        LOS north component.
    los_up : np.ndarray
        LOS up component.
    mask : np.ndarray
        Valid pixel mask.
    ref_point : tuple
        Reference pixel ``(row, col)``; the displacement is zeroed there
        before fitting and the offset is added back to the surface.
    window_size_x : int
        Window width
    window_size_y : int
        Window height
    wavelength_m : float
        Radar wavelength in meters, for unwrapping error correction.
    tropo_ref_file : Path, optional
        Per-epoch tropospheric correction for the reference date.
    tropo_sec_file : Path, optional
        Per-epoch tropospheric correction for the secondary date. The net
        ``sec - ref`` correction is removed before the fit.
    event_mask_file : Path, optional
        Per-epoch event mask GeoTIFF (1 = valid, 0 = event). The event region
        is filled from its neighbours before the fit, and the surface is
        still removed from the whole displacement. Without a mask, the
        residual-outlier options in `calibration_options` (if set) detect
        such regions automatically.
    fit_n_jobs : int, optional
        Parallel workers for the windowed fit, by default ``-1`` (all CPUs).
    fit_lock : threading.Lock, optional
        Held only during the windowed fit. `run` shares one lock across its
        per-file workers because a joblib ``Parallel`` fit nested inside
        another running one is silently throttled; serialising the fits
        lets each use every CPU while I/O still overlaps. ``None`` does not
        synchronise.
    recompute_gnss : bool, optional
        Ignore and overwrite the GNSS LOS caches, by default ``False``.
        `run` builds the frame-level caches once and passes ``False``.

    Returns
    -------
    Path or None
        Calibration surface path, or ``None`` if the file's dates cannot
        be parsed.

    """
    from ..gnss.reference import compute_gnss_los, compute_gnss_los_std
    from ..surface import estimate_calibration_surface
    from .utils import get_file_dates

    try:
        ref_date, sec_date = get_file_dates(disp_file)
    except Exception:
        logger.warning(f"Could not extract dates from {disp_file.name}")
        return None

    logger.debug(f"Processing {disp_file.name}")
    cal_opts = config.algorithm_parameters.calibration_options
    product_path = config.run_config.product_path_group.product_path

    disp = io_reader.read_netcdf(disp_file, variable="displacement").data

    corrections: list[np.ndarray] = []
    tropo_applied = tropo_ref_file is not None and tropo_sec_file is not None
    if tropo_applied:
        assert tropo_ref_file is not None
        assert tropo_sec_file is not None
        corrections.append(
            io_reader.read_geotiff(tropo_sec_file).data
            - io_reader.read_geotiff(tropo_ref_file).data
        )

    set_applied = False
    if cal_opts.apply_solid_earth_tide_correction:
        set_corr = io_reader.read_correction_layer(disp_file, "solid_earth_tide")
        if set_corr is None:
            logger.warning(
                f"No /corrections/solid_earth_tide in {disp_file.name}; "
                "calibrating without SET correction"
            )
        else:
            corrections.append(set_corr)
            set_applied = True

    gnss_kwargs = {
        "gnss_ref": gnss_manager,
        "los_east": los_east,
        "los_north": los_north,
        "los_up": los_up,
        "grid": disp_file,
        "cache_dir": product_path,
        "ref_date": ref_date,
        "sec_date": sec_date,
        "recompute": recompute_gnss,
    }
    # GNSS fields are in mm; the displacement is in m.
    gnss_los = compute_gnss_los(**gnss_kwargs) / 1000.0
    gnss_los_std = None
    if cal_opts.weight_fit_by_gnss_uncertainty:
        gnss_los_std = compute_gnss_los_std(**gnss_kwargs) / 1000.0

    event_mask = None
    if event_mask_file is not None:
        event_mask = io_reader.read_geotiff(event_mask_file).data.astype(bool)
        logger.debug(f"Event mask: {event_mask_file.name}")

    downsample_factor = config.grid_settings.downsample_factor
    weights = None
    if downsample_factor > 1 and config.grid_settings.downsample_weighted:
        weights = io_reader.read_netcdf(disp_file, variable="temporal_coherence").data
        logger.info("Downsampling weighted by temporal coherence")

    try:
        result = estimate_calibration_surface(
            disp,
            gnss_los,
            mask,
            ref_point,
            (window_size_x, window_size_y),
            corrections=corrections,
            gnss_los_std=gnss_los_std,
            event_mask=event_mask,
            options=cal_opts,
            wavelength_m=wavelength_m,
            downsample_factor=downsample_factor,
            downsample_method=config.grid_settings.downsample_method,
            downsample_weights=weights,
            n_jobs=fit_n_jobs,
            fit_lock=fit_lock,
            fit_surface=spatial_processor.fit_windowed_surface,
        )
    except ValueError as exc:
        exc.add_note(f"While calibrating {disp_file.name}")
        raise

    grid_type = config.grid_settings.grid_type
    ref_frame = config.grid_settings.reference_frame
    suffix = f"_calibration_surface_{grid_type}_{ref_frame.lower()}"
    if downsample_factor > 1:
        suffix += f"_downsample{downsample_factor}"
    if tropo_applied:
        suffix += "_tropo"
    if set_applied:
        suffix += "_set"
    if not cal_opts.unwrap_error_correction:
        suffix += "_nowrap"
    output_file = product_path / f"{disp_file.stem}{suffix}.tif"

    # The frame is recorded so a LOS decomposition can refuse to combine
    # geometries calibrated in different frames.
    description = f"LOS calibration surface ({grid_type} GNSS model, {ref_frame} frame)"
    if tropo_applied:
        description += " with tropospheric correction"
    if set_applied:
        description += " with solid Earth tide correction"
    if result.n_auto_masked_pixels:
        description += (
            f"; auto-masked {result.n_auto_masked_pixels} px "
            f"(MAD threshold={cal_opts.residual_outlier_mad_threshold})"
        )
    if result.n_region_masked_pixels:
        description += (
            f"; region-masked {result.n_region_masked_pixels} px "
            f"(region MAD threshold={cal_opts.residual_region_mad_threshold}, "
            f"min_pixels={cal_opts.residual_region_min_pixels})"
        )

    io_writer.write_geotiff(
        result.surface,
        output_file,
        reference_file=disp_file,
        nodata=np.nan,
        descriptions=[description],
    )

    logger.debug(f"Saved: {output_file.name}")
    return output_file


@dataclass
class CalibrationWorkflow:
    """Run the calibration step functions for one file or a whole stack.

    `run_single` calibrates one product, `run` the full batch.

    Attributes
    ----------
    config : WorkflowConfig
        Workflow configuration
    gnss_manager : GNSSReference
        GNSS reference data interface
    io_reader : RasterReader
        I/O reader for loading raster data
    io_writer : RasterWriter
        I/O writer for saving raster data
    spatial_processor : SpatialProcessor
        Spatial processor
    state : CalibrationState
        Workflow state tracking

    """

    config: WorkflowConfig
    gnss_manager: GNSSReference | None = field(default=None, init=False)
    io_reader: RasterReader = field(init=False)
    io_writer: RasterWriter = field(init=False)
    spatial_processor: SpatialProcessor = field(init=False)
    state: CalibrationState = field(default_factory=CalibrationState, init=False)

    def __post_init__(self):
        """Initialize workflow components and logging."""
        from ..io.read import RasterReader
        from ..io.write import RasterWriter
        from ..log_setup import configure_logging
        from ..spatial.processor import SpatialProcessor

        configure_logging(log_file=self.config.run_config.log_file)
        self.config.run_config.product_path_group.product_path.mkdir(
            parents=True, exist_ok=True
        )
        self.io_reader = RasterReader()
        self.io_writer = RasterWriter()
        self.spatial_processor = SpatialProcessor()

    def _setup(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[int, int]]:
        """Download GNSS stations, load LOS/mask, and select a reference point.

        Shared by `run_single` and `run`. Populates `self.gnss_manager`.

        Returns
        -------
        los_east, los_north, los_up, mask : np.ndarray
            See `load_los_and_mask`.
        ref_point : tuple of int
            See `find_reference_point`.

        """
        self.gnss_manager = setup_gnss_reference(
            io_reader=self.io_reader,
            input_files=self.config.input_options.input_files,
            product_path=self.config.run_config.product_path_group.product_path,
            reference_frame=self.config.grid_settings.reference_frame,
            grid_type=self.config.grid_settings.grid_type,
        )

        los_east, los_north, los_up, mask = load_los_and_mask(
            io_reader=self.io_reader,
            los_file=self.config.input_options.los_file,
            water_mask=self.config.input_options.water_mask,
            custom_mask=self.config.input_options.custom_mask,
        )

        ref_point = find_reference_point(
            input_files=self.config.input_options.input_files,
            product_path=self.config.run_config.product_path_group.product_path,
            mask=mask,
            reference_point=self.config.input_options.reference_point,
        )

        return los_east, los_north, los_up, mask, ref_point

    @with_scratch_temp_dir
    def run_single(
        self,
        disp_file: Path,
        tropo_ref_file: Path | None = None,
        tropo_sec_file: Path | None = None,
    ) -> CalibrationState:
        """Run the calibration workflow on a single displacement file.

        Performs the same setup as `run` (GNSS download, LOS/mask loading,
        reference point selection) but processes only the one specified file.

        Parameters
        ----------
        disp_file : Path
            Path to the NetCDF displacement file to calibrate.
        tropo_ref_file : Path, optional
            Per-epoch tropospheric correction GeoTIFF for the reference date.
        tropo_sec_file : Path, optional
            Per-epoch tropospheric correction GeoTIFF for the secondary date.

        Returns
        -------
        CalibrationState
            Workflow state with ``n_files_total = 1`` and, on success,
            one entry in ``output_files``.

        """
        from .utils import parse_window_size_meters

        logger.info("=" * 60)
        logger.info("Starting Venti Single-File Calibration")
        logger.info("=" * 60)
        logger.info(f"Input file : {disp_file}")

        if not self.config.grid_settings.apply_tropo_correction and (
            tropo_ref_file is not None or tropo_sec_file is not None
        ):
            logger.info("apply_tropo_correction=False; ignoring supplied tropo files")
            tropo_ref_file = None
            tropo_sec_file = None

        if tropo_ref_file is not None:
            logger.info(f"Tropo ref  : {tropo_ref_file}")
        if tropo_sec_file is not None:
            logger.info(f"Tropo sec  : {tropo_sec_file}")

        los_east, los_north, los_up, mask, ref_point = self._setup()
        gnss_manager = self.gnss_manager
        assert gnss_manager is not None, "gnss_manager not initialized"

        window_size_pixels = parse_window_size_meters(
            self.config.grid_settings.window_size_meters,
            self.config.grid_settings.posting_meters,
        )

        self.state.n_files_total = 1

        wavelength_m = self.config.input_options.wavelength_m
        logger.info("Radar wavelength: %.6f m", wavelength_m)

        event_mask_file = find_event_mask_file(
            self.config.input_options.event_mask_dir, disp_file
        )
        if event_mask_file is not None:
            logger.info(f"Event mask : {event_mask_file}")

        output_file = process_displacement_file(
            self.config,
            self.io_reader,
            self.io_writer,
            self.spatial_processor,
            gnss_manager,
            disp_file,
            los_east,
            los_north,
            los_up,
            mask,
            ref_point,
            window_size_pixels,
            window_size_pixels,
            wavelength_m=wavelength_m,
            tropo_ref_file=tropo_ref_file,
            tropo_sec_file=tropo_sec_file,
            event_mask_file=event_mask_file,
            recompute_gnss=self.config.grid_settings.recompute_gnss,
        )

        if output_file:
            self.state.output_files.append(output_file)
        else:
            self.state.n_files_failed += 1

        self.state.n_files_processed += 1

        logger.info("=" * 60)
        logger.info("Single-File Calibration Complete")
        logger.info("=" * 60)
        logger.info(f"Output: {output_file or 'FAILED'}")

        return self.state

    @with_scratch_temp_dir
    def run(self, max_files: int | None = None, n_workers: int = 1) -> CalibrationState:
        """Run the complete calibration workflow.

        Parameters
        ----------
        max_files : int or None, optional
            Processing only the first ``max_files`` displacement files.
            Default is None (process all files).
        n_workers : int, optional
            Files processed concurrently in threads, by default 1 (serial).
            Only I/O and pre-processing overlap: the windowed fits are
            serialised by a shared lock and each uses all CPUs (see
            `process_displacement_file`). 2-4 is a good range.

        Returns
        -------
        CalibrationState
            Final workflow state

        """
        logger.info("=" * 60)
        logger.info("Starting Venti Calibration Workflow")
        logger.info("=" * 60)

        los_east, los_north, los_up, mask, ref_point = self._setup()
        gnss_manager = self.gnss_manager
        assert gnss_manager is not None, "gnss_manager not initialized"

        disp_files = sorted(self.config.input_options.input_files.glob("*.nc"))
        if max_files is not None:
            disp_files = disp_files[:max_files]
        self.state.n_files_total = len(disp_files)

        logger.info(f"Processing {self.state.n_files_total} displacement files")

        from .utils import match_correction_to_displacement

        if (
            self.config.input_options.tropo_files is not None
            and self.config.grid_settings.apply_tropo_correction
        ):
            tropo_files = sorted(self.config.input_options.tropo_files.glob("*.tif"))
            matched_files = match_correction_to_displacement(tropo_files, disp_files)
            logger.info(f"Matched {len(matched_files)} tropospheric correction files")
        else:
            if self.config.input_options.tropo_files is not None:
                logger.info(
                    "tropo_files is configured but apply_tropo_correction=False; "
                    "skipping tropospheric correction"
                )
            matched_files = match_correction_to_displacement(None, disp_files)

        from .utils import parse_window_size_meters

        window_size_pixels = parse_window_size_meters(
            self.config.grid_settings.window_size_meters,
            self.config.grid_settings.posting_meters,
        )

        wavelength_m = self.config.input_options.wavelength_m
        logger.info("Radar wavelength: %.6f m", wavelength_m)

        fit_lock = threading.Lock() if n_workers > 1 else None

        if n_workers > 1:
            cpu_count = os.cpu_count() or 1
            logger.info(
                f"Parallel mode: {n_workers} file workers overlapping I/O; "
                f"surface fits serialised, each using all {cpu_count} CPUs "
                "per fit"
            )

        # 'constant' caches depend only on frame geometry: build them once so
        # epochs neither rebuild them nor race to write them in parallel.
        recompute_per_epoch = self.config.grid_settings.recompute_gnss
        if self.config.grid_settings.grid_type == "constant":
            from ..gnss.reference import compute_gnss_los, compute_gnss_los_std

            logger.info("Pre-computing GNSS LOS velocity cache...")
            gnss_cache_kwargs = {
                "gnss_ref": gnss_manager,
                "los_east": los_east,
                "los_north": los_north,
                "los_up": los_up,
                "grid": disp_files[0],
                "cache_dir": self.config.run_config.product_path_group.product_path,
                "recompute": self.config.grid_settings.recompute_gnss,
            }
            compute_gnss_los(**gnss_cache_kwargs)
            cal_opts = self.config.algorithm_parameters.calibration_options
            if cal_opts.weight_fit_by_gnss_uncertainty:
                compute_gnss_los_std(**gnss_cache_kwargs)
            recompute_per_epoch = False

        def _process_one(
            item: tuple[Path | None, Path | None, Path],
        ) -> Path | None:
            ref_tropo, sec_tropo, disp_file = item
            event_mask_file = find_event_mask_file(
                self.config.input_options.event_mask_dir, disp_file
            )
            return process_displacement_file(
                self.config,
                self.io_reader,
                self.io_writer,
                self.spatial_processor,
                gnss_manager,
                disp_file,
                los_east,
                los_north,
                los_up,
                mask,
                ref_point,
                window_size_pixels,
                window_size_pixels,
                wavelength_m=wavelength_m,
                tropo_ref_file=ref_tropo,
                tropo_sec_file=sec_tropo,
                event_mask_file=event_mask_file,
                fit_lock=fit_lock,
                recompute_gnss=recompute_per_epoch,
            )

        if n_workers > 1:
            from joblib import Parallel, delayed

            results: list[Path | None] = Parallel(n_jobs=n_workers, prefer="threads")(
                delayed(_process_one)(item) for item in matched_files
            )
        else:
            results = [
                _process_one(item) for item in tqdm(matched_files, desc="Calibrating")
            ]

        for output_file in results:
            if output_file:
                self.state.output_files.append(output_file)
            else:
                self.state.n_files_failed += 1
            self.state.n_files_processed += 1

        logger.info("=" * 60)
        logger.info("Calibration Workflow Complete!")
        logger.info("=" * 60)
        logger.info(f"Total files: {self.state.n_files_total}")
        logger.info(f"Processed: {self.state.n_files_processed}")
        logger.info(f"Failed: {self.state.n_files_failed}")
        logger.info(f"Success rate: {self.state.success_rate:.1f}%")
        logger.info(
            "Output directory:"
            f" {self.config.run_config.product_path_group.product_path}"
        )

        return self.state
