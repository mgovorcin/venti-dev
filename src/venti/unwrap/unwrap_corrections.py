"""Unwrap error correction module for interferometric phase unwrapping.

This module provides classes and functions to correct residual offsets in unwrapped
interferograms using watershed segmentation and regional median corrections.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, segmentation

# Optional imports for type hints
try:
    from rasterio.transform import Affine
except ImportError:
    Affine = None

# Set up logger
logger = logging.getLogger(__name__)


class UnwrapCorrector:
    """Class for correcting unwrapping errors in interferometric phase data.

    This class implements watershed-based segmentation to identify coherent
    regions and correct unwrap errors by detecting and removing integer
    multiples of the wavelength.

    Parameters
    ----------
    min_region_area : int, optional
        Minimum area (in pixels) for valid regions, by default 20
    wavelength : float, optional
        The wavelength value for unwrap correction, by default 0.0555m for Sentinel-1.
        For phase in radians, use 2π. For other satellites, use appropriate wavelength.

    Attributes
    ----------
    min_region_area : int
        Minimum area threshold for region filtering
    wavelength : float
        Wavelength used for unwrap correction
    labeled_regions_ : np.ndarray or None
        Labeled regions from watershed segmentation (set after correction)
    valid_labels_ : np.ndarray or None
        Valid region labels after area filtering (set after correction)
    medians_ : np.ndarray or None
        Median displacement values for each region (set after correction)
    unwrap_cycles_ : np.ndarray or None
        Number of wavelength cycles corrected for each region (set after correction)
    n_regions_ : int or None
        Number of valid regions identified (set after correction)

    Examples
    --------
    Correct Sentinel-1 displacement (meters) - default::

        corrector = UnwrapCorrector()
        corrected = corrector.correct(disp, mask)
        info = corrector.get_region_info()
        n_regions = info['n_regions']
        cycles = info['unwrap_cycles']
        print(f"Corrected {n_regions} regions with cycles: {cycles}")

    Correct phase data (radians)::

        corrector = UnwrapCorrector(wavelength=2*np.pi)
        corrected = corrector.correct(phase, mask)

    """

    def __init__(self, min_region_area: int = 20, wavelength: float = 0.0555):
        """Initialize the UnwrapCorrector.

        Parameters
        ----------
        min_region_area : int, optional
            Minimum area (in pixels) for valid regions, by default 20
        wavelength : float, optional
            The wavelength value for unwrap correction
            (default is 0.0555m for Sentinel-1).
            For phase in radians, use 2π. Common values:
            - Sentinel-1 (C-band): 0.0555 m
            - ALOS-2 (L-band): 0.236 m

        """
        self.min_region_area = min_region_area
        self.wavelength = wavelength
        self.labeled_regions_: np.ndarray | None = None
        self.valid_labels_: np.ndarray | None = None
        self.medians_: np.ndarray | None = None
        self.unwrap_cycles_: np.ndarray | None = None
        self.n_regions_: int | None = None

    def _prepare_displacement(
        self, input_disp: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Prepare displacement data for watershed segmentation.

        Parameters
        ----------
        input_disp : np.ndarray
            Input displacement field

        Returns
        -------
        scaled_disp : np.ndarray
            Scaled displacement array for watershed segmentation
        disp_mask : np.ndarray
            Mask of invalid displacement values

        """
        masked_disp = np.ma.masked_equal(input_disp, 0)
        scaled_disp = masked_disp - np.nanmin(masked_disp) + 1
        scaled_disp = np.nan_to_num(scaled_disp, nan=0)
        disp_mask = np.ma.masked_invalid(input_disp).mask

        return scaled_disp, disp_mask

    def _watershed_segmentation(
        self, scaled_disp: np.ndarray, mask: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Perform watershed segmentation to identify coherent regions.

        Parameters
        ----------
        scaled_disp : np.ndarray
            Scaled displacement array
        mask : np.ndarray
            Binary mask indicating valid data regions

        Returns
        -------
        labeled_regions : np.ndarray
            Labeled regions from watershed
        valid_labels : np.ndarray
            Valid region labels after area filtering

        """
        # Apply Sobel filter to create elevation map for watershed
        elevation_map = filters.sobel(mask)

        # Generate markers for watershed segmentation
        markers = np.zeros_like(scaled_disp, dtype=int)
        markers[scaled_disp < 1] = 1  # Marker for background
        markers[scaled_disp > 1] = 2  # Marker for objects

        # Perform watershed segmentation
        segmentation_result = segmentation.watershed(elevation_map, markers)
        segmentation_result = ndi.binary_fill_holes(segmentation_result - 1)
        labeled_regions, _ = ndi.label(segmentation_result)

        # Extract region properties and filter by area
        regions = measure.regionprops(labeled_regions)
        region_areas = np.array([r.area for r in regions])
        region_labels = np.array([r.label for r in regions])
        valid_labels = region_labels[region_areas > self.min_region_area]

        logger.info(f"Identified {len(valid_labels)} regions in watershed segmentation")

        return labeled_regions, valid_labels

    def _compute_regional_medians(
        self,
        input_disp: np.ndarray,
        labeled_regions: np.ndarray,
        valid_labels: np.ndarray,
    ) -> np.ndarray:
        """Compute median displacement for each region.

        Parameters
        ----------
        input_disp : np.ndarray
            Input displacement array
        labeled_regions : np.ndarray
            Labeled regions
        valid_labels : np.ndarray
            Valid region labels

        Returns
        -------
        np.ndarray
            Median displacement values for each region

        """
        medians = ndi.labeled_comprehension(
            input_disp,
            labeled_regions,
            index=valid_labels,
            func=np.nanmedian,
            out_dtype=np.float32,
            default=0,
        )
        return medians

    def _compute_unwrap_cycles(
        self, medians: np.ndarray, ref_label: int = 0
    ) -> np.ndarray:
        """Compute the number of wavelength cycles to correct for each region.

        This method calculates the offset of each region relative to the
        reference region and determines the nearest integer number of
        wavelength cycles that should be removed.

        Parameters
        ----------
        medians : np.ndarray
            Median displacement values for each region
        ref_label : int, optional
            Index of reference region, by default 0

        Returns
        -------
        np.ndarray
            Number of wavelength cycles to subtract from each region

        """
        ref_median = medians[ref_label]
        # Calculate offset relative to reference
        offsets = medians - ref_median

        # Compute the number of cycles (round to nearest integer)
        # Handle NaN values by replacing them with 0
        with np.errstate(invalid="ignore"):
            cycles = np.round(offsets / self.wavelength)
        cycles = np.nan_to_num(cycles, nan=0.0).astype(np.int32)

        logger.info(
            f"Unwrap correction: {np.count_nonzero(cycles)} of {cycles.size} "
            "regions shifted by whole wavelengths"
        )
        logger.debug(f"Unwrap cycles per region: {cycles}")
        logger.debug(f"Wavelength corrections (in units): {cycles * self.wavelength}")

        return cycles

    def _apply_corrections(
        self,
        input_disp: np.ndarray,
        labeled_regions: np.ndarray,
        valid_labels: np.ndarray,
        unwrap_cycles: np.ndarray,
        disp_mask: np.ndarray,
    ) -> np.ndarray:
        """Apply unwrap corrections by removing wavelength cycle offsets.

        Parameters
        ----------
        input_disp : np.ndarray
            Input displacement field
        labeled_regions : np.ndarray
            Labeled regions
        valid_labels : np.ndarray
            Valid region labels
        unwrap_cycles : np.ndarray
            Number of wavelength cycles to remove from each region
        disp_mask : np.ndarray
            Mask of invalid displacement values

        Returns
        -------
        np.ndarray
            Corrected displacement as a masked array; pixels outside the
            corrected regions are unchanged.

        """
        # One correction per label; pixels outside the corrected regions
        # (background, regions below min_region_area) get 0 and keep their
        # value. Writing only region pixels into an uninitialised array left
        # those pixels as leftover memory.
        corrections = np.zeros(labeled_regions.max() + 1, dtype=np.float32)
        corrections[valid_labels] = unwrap_cycles * self.wavelength
        disp_updated = (input_disp - corrections[labeled_regions]).astype(np.float32)
        return np.ma.masked_array(disp_updated, mask=disp_mask)

    def correct(
        self, input_disp: np.ndarray, mask: np.ndarray, ref_region: int = 0
    ) -> np.ndarray:
        """Correct unwrapping errors by removing wavelength cycle offsets.

        This method uses watershed segmentation to identify coherent regions
        based on the displacement field and corrects unwrap errors by detecting
        and removing integer multiples of the wavelength from each region.

        Parameters
        ----------
        input_disp : np.ndarray
            Input displacement field (unwrapped interferogram)
        mask : np.ndarray
            Binary mask indicating valid data regions
        ref_region : int, optional
            Index of reference region for correction, by default 0

        Returns
        -------
        np.ndarray
            Corrected displacement field as a masked array

        Notes
        -----
        The correction workflow:
        1. Prepare displacement data for segmentation
        2. Apply watershed segmentation to identify regions
        3. Calculate median displacement for each region
        4. Compute unwrap cycle offsets (n * wavelength)
        5. Remove cycle offsets from each region

        """
        # Prepare displacement data
        scaled_disp, disp_mask = self._prepare_displacement(input_disp)

        # Perform watershed segmentation
        labeled_regions, valid_labels = self._watershed_segmentation(scaled_disp, mask)

        if valid_labels.size == 0:
            logger.warning(
                "No valid regions found after watershed segmentation "
                f"(min_region_area={self.min_region_area}). "
                "Returning input displacement unchanged."
            )
            self.labeled_regions_ = labeled_regions
            self.valid_labels_ = valid_labels
            self.medians_ = np.array([], dtype=np.float32)
            self.unwrap_cycles_ = np.array([], dtype=np.int32)
            self.n_regions_ = 0
            return np.ma.masked_array(input_disp, mask=disp_mask)

        # Compute regional medians
        medians = self._compute_regional_medians(
            input_disp, labeled_regions, valid_labels
        )

        # Compute unwrap cycles for each region
        unwrap_cycles = self._compute_unwrap_cycles(medians, ref_label=ref_region)

        # Store results for inspection
        self.labeled_regions_ = labeled_regions
        self.valid_labels_ = valid_labels
        self.medians_ = medians
        self.unwrap_cycles_ = unwrap_cycles
        self.n_regions_ = len(valid_labels)

        # Apply corrections
        corrected_disp = self._apply_corrections(
            input_disp, labeled_regions, valid_labels, unwrap_cycles, disp_mask
        )

        return corrected_disp

    def get_region_info(self) -> dict | None:
        """Get information about identified regions.

        Returns
        -------
        dict or None
            Dictionary containing region information, or None if correction
            hasn't been run yet. Includes:
            - n_regions: Number of regions
            - valid_labels: Region labels
            - medians: Median displacement per region
            - unwrap_cycles: Number of wavelength cycles corrected per region
            - labeled_regions: Labeled region array

        """
        if self.n_regions_ is None:
            return None

        return {
            "n_regions": self.n_regions_,
            "valid_labels": self.valid_labels_,
            "medians": self.medians_,
            "unwrap_cycles": self.unwrap_cycles_,
            "labeled_regions": self.labeled_regions_,
            "wavelength": self.wavelength,
        }

    def save_geotiff(
        self,
        corrected_disp: np.ndarray,
        output_path: str | Path,
        reference_file: str | Path | None = None,
        transform: Affine | None = None,
        crs: str | None = None,
        nodata: float | None = None,
    ) -> None:
        """Save corrected displacement to a GeoTIFF file.

        Parameters
        ----------
        corrected_disp : np.ndarray
            Corrected displacement array to save
        output_path : str or Path
            Output path for the GeoTIFF file
        reference_file : str or Path, optional
            Reference raster file to copy georeferencing from
        transform : Affine, optional
            Affine transform for georeferencing (if reference_file not provided)
        crs : str, optional
            Coordinate reference system (if reference_file not provided)
        nodata : float, optional
            NoData value, by default np.nan

        Raises
        ------
        ValueError
            If neither reference_file nor (transform and crs) are provided
        ImportError
            If rasterio is not installed

        Examples
        --------
        Save with reference file for georeferencing::

            corrector = UnwrapCorrector()
            corrected = corrector.correct(disp, mask)
            corrector.save_geotiff(
                corrected, 'corrected.tif', reference_file='input.tif'
            )

        Or specify georeferencing manually::

            corrector.save_geotiff(
                corrected, 'corrected.tif', transform=transform, crs='EPSG:4326'
            )

        """
        from venti.io import write_geotiff

        write_geotiff(
            corrected_disp,
            output_path,
            transform=transform,
            crs=crs,
            nodata=nodata,
            reference_file=reference_file,
        )


def correct_region_offset(
    input_disp: np.ndarray | str | Path,
    mask: np.ndarray | None = None,
    wavelength: float = 0.0555,
    min_region_area: int = 20,
    output_file: str | Path | None = None,
) -> np.ndarray:
    """Correct unwrapping errors by removing wavelength cycle offsets.

    This is a convenience function that wraps the UnwrapCorrector class
    for simple use cases. It detects and removes integer multiples of the
    wavelength from different regions based on the displacement field.

    Parameters
    ----------
    input_disp : np.ndarray or str or Path
        Input displacement field (unwrapped interferogram) as array,
        or path to NetCDF file (.nc) to read. If NetCDF file, expects
        'displacement' and 'water_mask' variables.
    mask : np.ndarray, optional
        Binary mask indicating valid data regions. If None and input_disp
        is a NetCDF file, mask will be read from 'water_mask' variable.
        Required if input_disp is an array.
    wavelength : float, optional
        The wavelength value for unwrap correction, by default 0.0555m for Sentinel-1.
        Common values:
        - Sentinel-1 (C-band): 0.0555 m
        - ALOS-2 (L-band): 0.236 m
        - Phase data (radians): 2π (6.283)
    min_region_area : int, optional
        Minimum area (in pixels) for valid regions, by default 20
    output_file : str or Path, optional
        If provided, save the corrected displacement to this GeoTIFF file.
        Georeferencing is automatically extracted from input NetCDF file.

    Returns
    -------
    np.ndarray
        Corrected displacement field, NaN where invalid or masked

    Examples
    --------
    From NetCDF file (mask auto-loaded from 'water_mask' variable)::

        corrected = correct_region_offset(
            'displacement.nc',
            output_file='corrected.tif'
        )

    From array with explicit mask::

        corrected = correct_region_offset(disp, mask)

    From NetCDF with custom output::

        corrected = correct_region_offset(
            'displacement.nc',
            output_file='corrected.tif',
            wavelength=0.0555
        )

    For phase data (radians)::

        corrected = correct_region_offset(phase, mask, wavelength=2*np.pi)

    For ALOS-2 displacement (meters)::

        corrected = correct_region_offset(disp, mask, wavelength=0.236)

    Notes
    -----
    NetCDF files must contain:

    - 'displacement' variable: unwrapped displacement data
    - 'water_mask' variable: binary mask (1=valid, 0=invalid)
    - Coordinate information (x/y or lon/lat)

    For more control and access to intermediate results, use the
    UnwrapCorrector class directly::

        corrector = UnwrapCorrector(wavelength=0.0555, min_region_area=20)
        corrected = corrector.correct(disp, mask)
        corrector.save_geotiff(corrected, 'corrected.tif', reference_file='input.tif')
        region_info = corrector.get_region_info()
        print(f"Corrected {region_info['n_regions']} regions")
        print(f"Unwrap cycles: {region_info['unwrap_cycles']}")

    """
    # Check if input is a file path
    geo_info = None
    if isinstance(input_disp, str | Path):
        from venti.io import read_netcdf

        input_path = Path(input_disp)
        if input_path.suffix.lower() == ".nc":
            logger.info(f"Reading NetCDF file: {input_path}")
            disp_data, mask_data, geo_info = read_netcdf(
                input_path, variable="displacement", mask_variable="water_mask"
            )
            # Use mask from file if not provided
            if mask is None:
                mask = mask_data
        else:
            msg = (
                f"Unsupported file format: {input_path.suffix}. Use .nc for NetCDF"
                " files."
            )
            raise ValueError(msg)
    else:
        disp_data = input_disp
        if mask is None:
            msg = "mask parameter is required when input_disp is an array"
            raise ValueError(msg)

    # Apply mask to displacement data
    logger.info("Applying mask to displacement data...")

    # Create masked array - set invalid values to NaN
    disp_masked = np.where(mask.astype(bool), disp_data, np.nan)

    n_valid = np.sum(~np.isnan(disp_masked))
    n_total = disp_masked.size
    logger.info(f"Valid pixels: {n_valid}/{n_total} ({100 * n_valid / n_total:.1f}%)")

    # Run correction
    corrector = UnwrapCorrector(min_region_area=min_region_area, wavelength=wavelength)
    # Plain array with NaN for invalid pixels, like the rest of Venti's API.
    corrected = np.ma.filled(
        np.ma.asarray(corrector.correct(disp_masked, mask), dtype=float), np.nan
    )

    # Save to GeoTIFF if output file is specified
    if output_file is not None:
        from venti.io import write_geotiff

        if geo_info is None:
            msg = (
                "Cannot save to GeoTIFF without georeferencing information. Provide"
                " input_disp as a NetCDF file path to automatically extract"
                " georeferencing."
            )
            raise ValueError(msg)

        # Check that transform and crs were successfully extracted
        if "transform" not in geo_info or geo_info["transform"] is None:
            msg = (
                "Failed to extract transform from NetCDF file. "
                "Ensure the file has valid coordinate information (x/y or lon/lat)."
            )
            raise ValueError(msg)
        if "crs" not in geo_info or geo_info["crs"] is None:
            msg = (
                "Failed to extract CRS from NetCDF file. "
                "Ensure the file has valid coordinate reference system information."
            )
            raise ValueError(msg)

        write_geotiff(
            corrected,
            output_file,
            transform=geo_info["transform"],
            crs=geo_info["crs"],
            nodata=geo_info.get("nodata", np.nan),
        )

    return corrected
