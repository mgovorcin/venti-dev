# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""DEM generation for OPERA frames using GLO30 data.

Note: The DEM generated here might differ from the DEM in DISP-STATIC products,
which use GLO30 prepared specifically for NISAR processing.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pyproj
import tyro
import xarray as xr
from dem_stitcher import stitch_dem
from numpy.typing import NDArray
from opera_utils import get_frame_bbox

# Constants
DEFAULT_BUFFER_METERS = 10_000.0
DEM_SOURCE = "glo_30"
OUTPUT_CRS = "EPSG:4326"
DEM_NODATA = np.nan

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

# Suppress verbose logging from third-party libraries
logging.getLogger("rasterio").setLevel(logging.WARNING)
logging.getLogger("dem_stitcher").setLevel(logging.WARNING)


def _transform_bbox_with_buffer(
    bbox: tuple[float, float, float, float],
    src_epsg: int,
    buffer: float,
) -> tuple[float, float, float, float]:
    """Transform bounding box to WGS84 and apply buffer.

    Parameters
    ----------
    bbox : tuple[float, float, float, float]
        Bounding box as (left, bottom, right, top) in source CRS.
    src_epsg : int
        Source EPSG code.
    buffer : float
        Buffer distance in meters to apply before transformation.

    Returns
    -------
    tuple[float, float, float, float]
        Transformed bounding box as (min_lon, min_lat, max_lon, max_lat) in WGS84.

    """
    transformer = pyproj.Transformer.from_crs(
        f"EPSG:{src_epsg}", OUTPUT_CRS, always_xy=True
    )

    left, bottom, right, top = bbox
    min_lon, min_lat = transformer.transform(left - buffer, bottom - buffer)
    max_lon, max_lat = transformer.transform(right + buffer, top + buffer)

    return min_lon, min_lat, max_lon, max_lat


def _stitch_dem_data(
    bounds: tuple[float, float, float, float],
) -> tuple[NDArray[np.float32], dict]:
    """Stitch DEM data from GLO30 for given bounds.

    Parameters
    ----------
    bounds : tuple[float, float, float, float]
        Bounding box as (min_lon, min_lat, max_lon, max_lat) in WGS84.

    Returns
    -------
    dem : NDArray[np.float32]
        DEM elevation data array.
    metadata : dict
        Metadata including width, height, and transform.

    """
    dem, metadata = stitch_dem(
        bounds,
        dem_name=DEM_SOURCE,
        dst_ellipsoidal_height=True,
        dst_area_or_point="Point",
    )
    return dem, metadata


def _create_dem_dataarray(
    dem: NDArray[np.float32],
    metadata: dict,
) -> xr.DataArray:
    """Create an xarray DataArray from DEM data and metadata.

    Parameters
    ----------
    dem : NDArray[np.float32]
        DEM elevation data array.
    metadata : dict
        Metadata including width, height, and transform.

    Returns
    -------
    xr.DataArray
        DEM as an xarray DataArray with coordinates and CRS.

    """
    width = metadata["width"]
    height = metadata["height"]
    transform = metadata["transform"]

    x = np.linspace(transform.c, transform.c + transform.a * (width - 1), width)
    y = np.linspace(transform.f, transform.f + transform.e * (height - 1), height)

    dem_ds = xr.DataArray(
        dem,
        dims=("y", "x"),
        coords={"x": x, "y": y},
        attrs={"crs": OUTPUT_CRS, "nodata": DEM_NODATA},
    )
    dem_ds.rio.write_crs(OUTPUT_CRS, inplace=True)

    return dem_ds


def generate_frame_dem(
    frame_id: int,
    buffer: float = DEFAULT_BUFFER_METERS,
    output_dir: Path = Path("./dem"),
    use_disp_epsg: bool = False,
) -> Path:
    """Download and generate DEM for a given OPERA frame.

    Fetches GLO30 DEM data for the frame's bounding box with a buffer,
    saves it as a GeoTIFF in WGS84, and optionally reprojects to the
    frame's native EPSG coordinate system.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    buffer : float, optional
        Buffer distance in meters around the frame extent.
        Default is 10 km (10,000 meters).
    output_dir : Path, optional
        Directory to save the resulting DEM GeoTIFF.
        Default is "./dem".
    use_disp_epsg : bool, optional
        If True, also save DEM reprojected to the frame's native EPSG.
        Default is False.

    Returns
    -------
    Path
        Path to the generated DEM file in WGS84.

    Raises
    ------
    ValueError
        If frame bbox cannot be retrieved.

    Examples
    --------
    >>> dem_path = generate_frame_dem(frame_id=8887, buffer=10000.0)
    >>> print(dem_path)
    ./dem/dem_frame_8887.tif

    """
    logger.info(f"Generating GLO30 DEM for frame {frame_id}")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"dem_frame_{frame_id}.tif"

    if output_path.exists():
        logger.info(f"DEM already exists: {output_path}")
        return output_path

    epsg, bbox = get_frame_bbox(frame_id)
    min_lon, min_lat, max_lon, max_lat = _transform_bbox_with_buffer(
        (bbox.left, bbox.bottom, bbox.right, bbox.top),
        epsg,
        buffer,
    )

    logger.info(
        "Fetching DEM for bounds: "
        f"[{min_lon:.4f}, {min_lat:.4f}, {max_lon:.4f}, {max_lat:.4f}]"
    )

    dem, metadata = _stitch_dem_data((min_lon, min_lat, max_lon, max_lat))
    dem_ds = _create_dem_dataarray(dem, metadata)

    dem_ds.rio.to_raster(output_path)
    logger.info(f"DEM saved to {output_path}")

    if use_disp_epsg:
        utm_path = output_dir / f"dem_frame_{frame_id}_epsg{epsg}.tif"
        dem_utm = dem_ds.rio.reproject(f"EPSG:{epsg}")
        dem_utm.rio.to_raster(utm_path)
        logger.info(f"DEM in EPSG:{epsg} saved to {utm_path}")

    return output_path


def main(
    frame_id: int,
    buffer: float = DEFAULT_BUFFER_METERS,
    output_dir: Path = Path("./dem"),
    use_disp_epsg: bool = False,
) -> None:
    """CLI entry point for DEM generation.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    buffer : float, optional
        Buffer distance in meters around frame extent. Default is 10 km.
    output_dir : Path, optional
        Output directory for DEM files. Default is "./dem".
    use_disp_epsg : bool, optional
        Whether to also generate DEM in frame's native EPSG. Default is False.

    """
    dem_path = generate_frame_dem(frame_id, buffer, output_dir, use_disp_epsg)
    logger.info(f"DEM generation complete: {dem_path}")


if __name__ == "__main__":
    tyro.cli(main)
