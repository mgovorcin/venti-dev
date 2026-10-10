# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""DISP-S1 CLI for generating LOS ENU and incidence angle rasters.

This module provides functionality to download CSLC-STATIC products from ASF,
stitch geometry layers, and generate line-of-sight (LOS) unit vectors in
East-North-Up (ENU) coordinates and incidence angle rasters for OPERA frames.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import asf_search as asf
import numpy as np
import rasterio as rio
import tyro
from numpy.typing import NDArray
from opera_utils import get_frame_to_burst_mapping
from opera_utils.geometry import stitch_geometry_layers
from osgeo import gdal

# Constants
PROCESSING_LEVEL = "CSLC-STATIC"
DEFAULT_DOWNLOAD_PROCESSES = 5
TEMP_DIR_NAME = "tmp"

# GeoTIFF output settings
GEOTIFF_DRIVER = "GTiff"
GEOTIFF_DTYPE = "float32"
COMPRESSION_LEVEL = 4
TILE_SIZE = 128

# Band indices (1-indexed for rasterio)
BAND_EAST = 1
BAND_NORTH = 2
BAND_UP = 3

# LOS component names for VRT generation
LOS_COMPONENTS = ["los_east", "los_north", "los_up"]

logger = logging.getLogger("disp_los_cli")
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logging.getLogger("asf_search").setLevel(logging.WARNING)


@contextmanager
def _temporary_directory(base_dir: Path) -> Iterator[Path]:
    """Create and cleanup a temporary directory.

    Parameters
    ----------
    base_dir : Path
        Base directory where temporary directory will be created.

    Yields
    ------
    Path
        Path to the temporary directory.

    """
    temp_dir = base_dir / TEMP_DIR_NAME
    temp_dir.mkdir(parents=True, exist_ok=True)
    try:
        yield temp_dir
    finally:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
            logger.info(f"Cleaned up temporary directory: {temp_dir}")


def _download_cslc_static_products(frame_id: int, output_dir: Path) -> list[Path]:
    """Download CSLC-STATIC products for a given frame.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Directory where products will be downloaded.

    Returns
    -------
    list[Path]
        Paths to downloaded CSLC-STATIC files.

    Raises
    ------
    ValueError
        If no burst IDs found for the frame or no ASF results returned.

    """
    burst_ids = get_frame_to_burst_mapping(frame_id).get("burst_ids", [])
    if not burst_ids:
        msg = f"No burst IDs found for frame {frame_id}"
        raise ValueError(msg)

    results = asf.search(operaBurstID=list(burst_ids), processingLevel=PROCESSING_LEVEL)
    if not results:
        msg = f"No {PROCESSING_LEVEL} products found for frame {frame_id}"
        raise ValueError(msg)

    logger.info(f"Downloading {len(results)} {PROCESSING_LEVEL} files to {output_dir}")
    results.download(path=output_dir, processes=DEFAULT_DOWNLOAD_PROCESSES)

    return [Path(output_dir, r.properties["fileName"]) for r in results]


def _load_los_components(
    temp_dir: Path,
) -> tuple[NDArray[np.float32], NDArray[np.float32], rio.Affine, rio.CRS]:
    """Load LOS east and north components from stitched geometry layers.

    Parameters
    ----------
    temp_dir : Path
        Directory containing stitched LOS rasters.

    Returns
    -------
    los_east : NDArray[np.float32]
        LOS east component.
    los_north : NDArray[np.float32]
        LOS north component.
    transform : rio.Affine
        Geospatial transform.
    crs : rio.CRS
        Coordinate reference system.

    Raises
    ------
    FileNotFoundError
        If LOS component files are missing after stitching.

    """
    los_east_path = temp_dir / "los_east.tif"
    los_north_path = temp_dir / "los_north.tif"

    if not (los_east_path.exists() and los_north_path.exists()):
        msg = (
            "Missing LOS east or north files after stitching. "
            f"Expected: {los_east_path}, {los_north_path}"
        )
        raise FileNotFoundError(msg)

    with rio.open(los_east_path) as east_ds, rio.open(los_north_path) as north_ds:
        los_east = east_ds.read(1).astype(np.float32)
        los_north = north_ds.read(1).astype(np.float32)
        transform = east_ds.transform
        crs = east_ds.crs

    return los_east, los_north, transform, crs


def _compute_los_up(
    los_east: NDArray[np.float32], los_north: NDArray[np.float32]
) -> NDArray[np.float32]:
    """Compute LOS up component from east and north components.

    Parameters
    ----------
    los_east : NDArray[np.float32]
        LOS east component.
    los_north : NDArray[np.float32]
        LOS north component.

    Returns
    -------
    NDArray[np.float32]
        LOS up component computed as sqrt(1 - east^2 - north^2).

    """
    return np.sqrt(1 - los_east**2 - los_north**2).astype(np.float32)


def _create_geotiff_profile(
    height: int,
    width: int,
    crs: rio.CRS,
    transform: rio.Affine,
    count: int = 3,
) -> dict:
    """Create a rasterio profile for GeoTIFF output.

    Parameters
    ----------
    height : int
        Raster height in pixels.
    width : int
        Raster width in pixels.
    crs : rio.CRS
        Coordinate reference system.
    transform : rio.Affine
        Geospatial transform.
    count : int, optional
        Number of bands, by default 3.

    Returns
    -------
    dict
        Rasterio profile dictionary.

    """
    return {
        "driver": GEOTIFF_DRIVER,
        "height": height,
        "width": width,
        "count": count,
        "dtype": GEOTIFF_DTYPE,
        "crs": crs,
        "transform": transform,
        "compress": "deflate",
        "zlevel": COMPRESSION_LEVEL,
        "tiled": True,
        "blockxsize": TILE_SIZE,
        "blockysize": TILE_SIZE,
        "predictor": 2,
        "interleave": "pixel",
    }


def _write_los_enu_raster(
    output_path: Path,
    los_east: NDArray[np.float32],
    los_north: NDArray[np.float32],
    los_up: NDArray[np.float32],
    crs: rio.CRS,
    transform: rio.Affine,
) -> None:
    """Write LOS ENU components to a 3-band GeoTIFF.

    Parameters
    ----------
    output_path : Path
        Output file path.
    los_east : NDArray[np.float32]
        LOS east component.
    los_north : NDArray[np.float32]
        LOS north component.
    los_up : NDArray[np.float32]
        LOS up component.
    crs : rio.CRS
        Coordinate reference system.
    transform : rio.Affine
        Geospatial transform.

    """
    height, width = los_east.shape
    profile = _create_geotiff_profile(height, width, crs, transform, count=3)

    with rio.open(output_path, "w", **profile) as dst:
        dst.write(los_east, BAND_EAST)
        dst.write(los_north, BAND_NORTH)
        dst.write(los_up, BAND_UP)
        dst.set_band_description(BAND_EAST, "LOS East")
        dst.set_band_description(BAND_NORTH, "LOS North")
        dst.set_band_description(BAND_UP, "LOS Up")

    logger.info(f"LOS ENU raster saved to {output_path}")


def _create_component_vrts(output_path: Path) -> None:
    """Create VRT files for each LOS component band.

    Parameters
    ----------
    output_path : Path
        Path to the multi-band LOS ENU raster.

    """
    for band, name in enumerate(LOS_COMPONENTS, start=1):
        vrt_path = output_path.parent / f"{name}.vrt"
        gdal.BuildVRT(str(vrt_path), [str(output_path)], bandList=[band])
        logger.info(f"VRT for {name} saved to {vrt_path}")


def generate_los_enu_raster(frame_id: int, output_dir: Path) -> Path:
    """Generate LOS ENU raster for a given OPERA frame.

    Downloads CSLC-STATIC products, stitches geometry layers, computes the
    LOS up component, and writes a 3-band GeoTIFF with ENU components.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Directory where output raster will be saved.

    Returns
    -------
    Path
        Path to the generated LOS ENU raster.

    Raises
    ------
    ValueError
        If no burst IDs or ASF results found for the frame.
    FileNotFoundError
        If stitched LOS component files are missing.

    Examples
    --------
    >>> los_path = generate_los_enu_raster(frame_id=123, output_dir=Path("./los"))
    >>> print(los_path)
    ./los/los_enu_frame_123.tif

    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"los_enu_frame_{frame_id}.tif"

    if output_path.exists():
        logger.info(f"LOS raster already exists: {output_path}")
        return output_path

    with _temporary_directory(output_dir) as temp_dir:
        static_files = _download_cslc_static_products(frame_id, temp_dir)
        stitch_geometry_layers(static_files, output_dir=temp_dir)

        los_east, los_north, transform, crs = _load_los_components(temp_dir)
        los_up = _compute_los_up(los_east, los_north)

        _write_los_enu_raster(output_path, los_east, los_north, los_up, crs, transform)
        _create_component_vrts(output_path)

    return output_path


def generate_incidence_angle_raster(los_enu_path: Path, output_dir: Path) -> Path:
    """Generate incidence angle raster from LOS ENU raster.

    Computes the incidence angle as arccos(los_up) and saves as a single-band
    GeoTIFF in degrees.

    Parameters
    ----------
    los_enu_path : Path
        Path to the LOS ENU raster.
    output_dir : Path
        Directory where output raster will be saved.

    Returns
    -------
    Path
        Path to the generated incidence angle raster.

    Raises
    ------
    FileNotFoundError
        If the input LOS ENU file does not exist.
    ValueError
        If frame ID cannot be extracted from filename.

    Examples
    --------
    >>> inc_path = generate_incidence_angle_raster(
    ...     los_enu_path=Path("los_enu_frame_123.tif"),
    ...     output_dir=Path("./los")
    ... )
    >>> print(inc_path)
    ./los/incidence_angle_frame_123.tif

    """
    if not los_enu_path.exists():
        msg = f"LOS ENU file does not exist: {los_enu_path}"
        raise FileNotFoundError(msg)

    frame_id = los_enu_path.stem.split("_")[-1]
    if not frame_id.isdigit():
        msg = f"Cannot extract frame ID from filename: {los_enu_path.name}"
        raise ValueError(msg)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"incidence_angle_frame_{frame_id}.tif"

    if output_path.exists():
        logger.info(f"Incidence angle raster already exists: {output_path}")
        return output_path

    with rio.open(los_enu_path) as los_ds:
        los_up = los_ds.read(BAND_UP).astype(np.float32)
        transform = los_ds.transform
        crs = los_ds.crs
        height, width = los_up.shape

    # Correct NaN masking - cannot use == with NaN
    mask = np.isnan(los_up)
    los_up_masked = np.ma.masked_array(los_up, mask=mask)
    incidence_angle = np.rad2deg(np.arccos(los_up_masked)).filled(0).astype(np.float32)

    profile = _create_geotiff_profile(height, width, crs, transform, count=1)

    with rio.open(output_path, "w", **profile) as dst:
        dst.write(incidence_angle, 1)
        dst.set_band_description(1, "Incidence Angle (degrees)")

    logger.info(f"Incidence angle raster saved to {output_path}")
    return output_path


@dataclass
class GenerateLOS:
    """CLI for generating LOS ENU and incidence angle rasters.

    Attributes
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Directory where output rasters will be saved.
    inc_out : bool
        Whether to generate incidence angle raster in addition to LOS ENU.

    """

    frame_id: int
    output_dir: Path = Path("./los")
    inc_out: bool = True

    def __call__(self) -> None:
        """Execute LOS raster generation."""
        los_path = generate_los_enu_raster(self.frame_id, self.output_dir)
        if self.inc_out:
            generate_incidence_angle_raster(los_path, self.output_dir)


def main() -> None:
    """Run CLI."""
    tyro.cli(GenerateLOS)()


if __name__ == "__main__":
    main()
