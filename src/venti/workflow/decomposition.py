"""Object-oriented decomposition workflow for InSAR displacement products.

This module provides a high-level DecompositionWorkflow class that orchestrates
the LOS-to-ENU (Line-of-Sight to East-North-Up) decomposition process.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from tqdm import tqdm

from .config import WorkflowConfig
from .utils import with_scratch_temp_dir

if TYPE_CHECKING:
    from ..io.read import RasterReader
    from ..io.write import RasterWriter
    from ..spatial.processor import SpatialProcessor

logger = logging.getLogger(__name__)


@dataclass
class DecompositionState:
    """State tracking for decomposition workflow.

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


@dataclass
class DecompositionWorkflow:
    """High-level workflow manager for InSAR LOS-to-ENU decomposition.

    This class orchestrates the decomposition process, projecting Line-of-Sight
    displacement measurements from multiple viewing geometries (ascending/descending)
    into East-North-Up displacement components.

    Attributes
    ----------
    config : WorkflowConfig
        Workflow configuration
    io_reader : RasterReader
        I/O reader for loading raster data
    io_writer : RasterWriter
        I/O writer for saving raster data
    spatial_processor : SpatialProcessor
        Spatial processor for decomposition operations
    state : DecompositionState
        Workflow state tracking

    """

    config: WorkflowConfig
    io_reader: RasterReader | None = field(default=None, init=False)
    io_writer: RasterWriter | None = field(default=None, init=False)
    spatial_processor: SpatialProcessor | None = field(default=None, init=False)
    state: DecompositionState = field(default_factory=DecompositionState, init=False)

    def __post_init__(self):
        """Initialize workflow components and logging."""
        from ..io.read import RasterReader
        from ..io.write import RasterWriter
        from ..log_setup import configure_logging
        from ..spatial.processor import SpatialProcessor

        configure_logging(log_file=self.config.run_config.log_file)

        # Create output directory
        self.config.input_options.work_directory.mkdir(parents=True, exist_ok=True)

        # Initialize components
        self.io_reader = RasterReader()
        self.io_writer = RasterWriter()
        self.spatial_processor = SpatialProcessor()

        # Store worker settings for use in processing methods
        self.gpu_enabled = self.config.worker_settings.gpu_enabled
        self.block_shape = self.config.worker_settings.block_shape
        self.threads_per_worker = self.config.worker_settings.threads_per_worker

        logger.info("Decomposition workflow initialized")

    def load_los_vectors(self) -> dict[str, np.ndarray]:
        """Load LOS vectors for all geometries.

        Returns
        -------
        dict
            Dictionary with 'ascending' and 'descending' LOS vectors
            Each contains (los_east, los_north, los_up) arrays

        """
        assert self.io_reader is not None, "io_reader not initialized"

        logger.info("Loading LOS unit vectors...")

        los_vectors = {}

        # Load LOS for primary geometry (from config)
        los_data = self.io_reader.read_geotiff(self.config.input_options.los_file)
        if los_data.data.ndim == 3:
            los_east = los_data.data[0]
            los_north = los_data.data[1]
            los_up = los_data.data[2]
        else:
            msg = "LOS file must be 3-band (east, north, up)"
            raise ValueError(msg)

        # TODO: Determine geometry type (ascending/descending) from metadata
        # For now, placeholder assumes ascending
        los_vectors["primary"] = (los_east, los_north, los_up)

        logger.info("LOS vectors loaded")
        return los_vectors

    def load_displacement_pairs(self) -> list[tuple[Path, Path | None]]:
        """Load and match displacement files from multiple geometries.

        Returns
        -------
        list of tuple
            List of (ascending_file, descending_file) pairs matched by date

        Notes
        -----
        This is a placeholder. Implementation should:
        1. Find displacement files from ascending and descending directories
        2. Match files by acquisition date
        3. Return paired files for decomposition

        """
        logger.info("Loading displacement file pairs...")

        # Placeholder: Load from single directory
        disp_files = sorted(self.config.input_options.input_files.glob("*.nc"))
        self.state.n_files_total = len(disp_files)

        # TODO: Implement proper pairing logic for ascending/descending
        # For now, create dummy pairs
        pairs: list[tuple[Path, Path | None]] = [(f, None) for f in disp_files]

        logger.info(f"Found {len(pairs)} displacement file pairs")
        return pairs

    def decompose_to_enu(
        self,
        los_asc: np.ndarray,
        _los_desc: np.ndarray,
        los_vectors_asc: tuple[np.ndarray, np.ndarray, np.ndarray],
        los_vectors_desc: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Decompose LOS displacements to East-North-Up components.

        Parameters
        ----------
        los_asc : np.ndarray
            LOS displacement from ascending geometry
        los_desc : np.ndarray
            LOS displacement from descending geometry
        los_vectors_asc : tuple of np.ndarray
            (east, north, up) unit vectors for ascending
        los_vectors_desc : tuple of np.ndarray
            (east, north, up) unit vectors for descending

        Returns
        -------
        tuple of np.ndarray
            (east_displacement, north_displacement, up_displacement)

        Notes
        -----
        This is a placeholder. Full implementation should:
        1. Set up least-squares inversion: d = G * m
           where d = [los_asc, los_desc], m = [east, north, up]
        2. Apply proper weighting based on geometry and uncertainties
        3. Handle rank-deficient cases (e.g., only ascending or descending)
        4. Compute uncertainties for output components

        """
        logger.debug("Performing LOS-to-ENU decomposition...")

        # Placeholder: Simple 2D decomposition (East-Up only)
        # Assumes purely vertical and horizontal motion
        # This is NOT the proper implementation!

        _e_asc, _n_asc, _u_asc = los_vectors_asc
        _e_desc, _n_desc, _u_desc = los_vectors_desc

        # TODO: Implement proper least-squares inversion
        # For now, return placeholder arrays
        shape = los_asc.shape
        east_disp = np.zeros(shape)
        north_disp = np.zeros(shape)
        up_disp = np.zeros(shape)

        logger.warning("Using placeholder decomposition - implement proper inversion!")

        return east_disp, north_disp, up_disp

    def process_displacement_pair(
        self,
        asc_file: Path,
        _desc_file: Path | None,
        _los_vectors: dict[str, tuple],
    ) -> tuple[Path, Path, Path] | None:
        """Process a pair of displacement files to produce ENU components.

        Parameters
        ----------
        asc_file : Path
            Ascending geometry displacement file
        desc_file : Path or None
            Descending geometry displacement file
        los_vectors : dict
            Dictionary of LOS vectors for each geometry

        Returns
        -------
        tuple of Path or None
            (east_file, north_file, up_file) output paths if successful

        Notes
        -----
        This is a placeholder. Full implementation should:
        1. Read both ascending and descending displacement files
        2. Apply any necessary corrections (atmosphere, reference point, etc.)
        3. Perform decomposition to ENU
        4. Write three output files (East, North, Up components)
        5. Include proper error handling and quality metrics

        """
        assert self.io_reader is not None, "io_reader not initialized"
        assert self.io_writer is not None, "io_writer not initialized"

        logger.debug(f"Processing {asc_file.name}")

        try:
            # Read displacement data
            netcdf_data = self.io_reader.read_netcdf(asc_file, variable="displacement")
            disp = netcdf_data.data

            # Placeholder: Create dummy ENU outputs
            east_disp = np.zeros_like(disp)
            north_disp = np.zeros_like(disp)
            up_disp = disp  # Placeholder: assume all displacement is vertical

            # Build output filenames
            output_dir = self.config.input_options.work_directory
            base_name = asc_file.stem

            east_file = output_dir / f"{base_name}_east.tif"
            north_file = output_dir / f"{base_name}_north.tif"
            up_file = output_dir / f"{base_name}_up.tif"

            # Save components
            self.io_writer.write_geotiff(
                east_disp,
                east_file,
                reference_file=asc_file,
                nodata=np.nan,
                descriptions=["East displacement component"],
            )

            self.io_writer.write_geotiff(
                north_disp,
                north_file,
                reference_file=asc_file,
                nodata=np.nan,
                descriptions=["North displacement component"],
            )

            self.io_writer.write_geotiff(
                up_disp,
                up_file,
                reference_file=asc_file,
                nodata=np.nan,
                descriptions=["Up displacement component"],
            )

        except Exception:
            logger.exception(f"Failed to process {asc_file.name}")
            return None
        else:
            logger.debug(f"Saved: {east_file.name}, {north_file.name}, {up_file.name}")
            return east_file, north_file, up_file

    @with_scratch_temp_dir
    def run(self) -> DecompositionState:
        """Run the complete decomposition workflow.

        Returns
        -------
        DecompositionState
            Final workflow state

        Notes
        -----
        Current implementation is a PLACEHOLDER. Full workflow should:
        1. Load LOS vectors for all available geometries
        2. Load and pair displacement files from multiple geometries
        3. For each pair:
           - Read displacement data
           - Perform LOS-to-ENU decomposition
           - Save East, North, Up components
        4. Generate quality metrics and uncertainty estimates
        5. Optionally create visualization products

        """
        logger.info("=" * 60)
        logger.info("Starting Venti Decomposition Workflow")
        logger.info("=" * 60)
        logger.warning(
            "PLACEHOLDER IMPLEMENTATION - Full decomposition not yet implemented!"
        )

        # Load LOS vectors
        los_vectors = self.load_los_vectors()

        # Load displacement file pairs
        disp_pairs = self.load_displacement_pairs()

        # Process each pair
        for asc_file, desc_file in tqdm(disp_pairs, desc="Decomposing"):
            output_files = self.process_displacement_pair(
                asc_file,
                desc_file,
                los_vectors,
            )

            if output_files:
                self.state.output_files.extend(output_files)
            else:
                self.state.n_files_failed += 1

            self.state.n_files_processed += 1

        # Summary
        logger.info("=" * 60)
        logger.info("Decomposition Workflow Complete!")
        logger.info("=" * 60)
        logger.info(f"Total pairs: {self.state.n_files_total}")
        logger.info(f"Processed: {self.state.n_files_processed}")
        logger.info(f"Failed: {self.state.n_files_failed}")
        logger.info(f"Success rate: {self.state.success_rate:.1f}%")
        logger.info(f"Output directory: {self.config.input_options.work_directory}")

        return self.state
