# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""DISP Frame Model Generator.

This script generates GIA and plate motion model grids for specific DISP
frames. It creates raster files showing horizontal and vertical velocity
components in mm/yr for the specified tectonic plate in ITRF 2014 or 2020
realization and for specific GIA model either Caron et al 2018 or ICE6D
(Peltier et al 2018).

Example usage:
    python get_model_rates.py 1234 --output results --plate EURA --itrf 2020 --gia ICE6
    python get_model_rates.py 5678 --output /path/to/output --plate PCFC
"""

import logging
from contextlib import suppress
from pathlib import Path
from typing import Any, Literal, NamedTuple

import geopandas as gpd
import numpy as np
import opera_utils
import rasterio
import tyro
from geepers.euler import plate_velocity_enu
from pyproj import Transformer
from rasterio.transform import Affine, from_bounds
from rasterio.warp import Resampling, reproject
from tyro import conf

from venti.models import load_gia

# Constants
DEFAULT_CELL_SIZE = 30
DEFAULT_GRID_POSTING = 10000
DEFAULT_PMM_GRID_POSTING = 10000
SUPPORTED_ITRF_YEARS = [2014, 2020]
SUPPORTED_GIA_MODELS = ["CARON", "ICE6"]

logger = logging.getLogger("venti")


class FrameUTM(NamedTuple):
    """Represents a UTM frame with coordinate system.

    Attributes
    ----------
    crs : object
        Coordinate reference system of the UTM frame.
    bounds : tuple
        Bounding box coordinates of the frame (minx, miny, maxx, maxy).
    transform : Affine
        Affine transform for rasterization.
    width : int
        Width of the raster grid in pixels.
    height : int
        Height of the raster grid in pixels.

    """

    crs: Any
    bounds: Any
    transform: Any
    width: int
    height: int


def get_frame_utm_transform(
    frame_gdf: gpd.GeoDataFrame, cell_size: int = DEFAULT_CELL_SIZE
) -> FrameUTM:
    """Compute UTM transform and grid dimensions for a GeoDataFrame frame.

    Parameters
    ----------
    frame_gdf : GeoDataFrame
        GeoDataFrame containing the frame geometry.
    cell_size : int, optional
        Grid cell size in meters (default is 30).

    Returns
    -------
    FrameUTM
        NamedTuple containing CRS, bounds, transform, width, and height.

    """
    utm_crs = frame_gdf.estimate_utm_crs()
    selected_frame_utm = frame_gdf.to_crs(utm_crs)
    minx, miny, maxx, maxy = selected_frame_utm.total_bounds

    xs = np.arange(minx, maxx, cell_size)
    ys = np.arange(maxy, miny, -cell_size)

    dst_transform = Affine(cell_size, 0, minx, 0, -cell_size, maxy)
    target_height, target_width = len(ys), len(xs)

    return FrameUTM(
        crs=utm_crs,
        bounds=selected_frame_utm.total_bounds,
        transform=dst_transform,
        width=target_width,
        height=target_height,
    )


def make_point_grid(
    gdf: gpd.GeoDataFrame,
    cell_size: float = DEFAULT_GRID_POSTING,
    target_crs: str = "EPSG:4326",
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a regular point grid within a polygon GeoDataFrame.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input polygon(s), assumed to be in EPSG:4326.
    cell_size : float, optional
        Spacing between points in meters (default is 50,000 = 50 km).
    target_crs : str, optional
        CRS for output grid (default is "EPSG:4326").

    Returns
    -------
    x : ndarray
        X coordinates of grid points (longitude if EPSG:4326).
    y : ndarray
        Y coordinates of grid points (latitude if EPSG:4326).

    Raises
    ------
    ValueError
        If input GeoDataFrame is empty or invalid.

    """
    if gdf.empty:
        msg = "Input GeoDataFrame is empty"
        raise ValueError(msg)

    if gdf.crs is None:
        logger.warning("No CRS defined for input GeoDataFrame, assuming EPSG:4326")
        gdf = gdf.set_crs("EPSG:4326")

    try:
        # Project polygon to UTM for accurate meter-based calculations
        utm_crs = gdf.estimate_utm_crs()
        gdf_proj = gdf.to_crs(utm_crs)
        logger.info(f"Projected to {utm_crs} for grid generation")

        # Get bounds
        minx, miny, maxx, maxy = gdf_proj.total_bounds
        logger.info(f"Bounds: {minx:.0f}, {miny:.0f}, {maxx:.0f}, {maxy:.0f}")

        # Calculate grid dimensions
        width = max(1, int(np.rint((maxx - minx) / cell_size)))
        length = max(1, int(np.rint((maxy - miny) / cell_size)))
        logger.info(f"Grid dimensions: {width} x {length} points")

        # Generate grid coordinates
        gridy, gridx = np.mgrid[
            maxy : miny : complex(length), minx : maxx : complex(width)  # type: ignore[misc]
        ]

        # Transform back to target CRS
        transformer = Transformer.from_crs(utm_crs, target_crs, always_xy=True)
        x, y = transformer.transform(gridx, gridy)

        logger.info(f"Generated grid with {x.size} points")
    except Exception:
        logger.exception("Error generating grid")
        raise
    else:
        return x, y


def _validate_itrf_year(date: int) -> None:
    """Validate the ITRF reference year.

    Parameters
    ----------
    date : int
        ITRF reference year.

    Raises
    ------
    ValueError
        If an unsupported ITRF year is specified.

    """
    if date not in SUPPORTED_ITRF_YEARS:
        msg = f"ITRF year must be one of {SUPPORTED_ITRF_YEARS}, got {date}"
        raise ValueError(msg)


def get_frame_pmm(
    frame_gdf: gpd.GeoDataFrame,
    plate: str = "NOAM",
    date: int = 2014,
    grid_posting: float = DEFAULT_PMM_GRID_POSTING,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Generate plate motion model (PMM) velocities for a given frame.

    Parameters
    ----------
    frame_gdf : GeoDataFrame
        GeoDataFrame representing the area of interest.
    plate : str, optional
        Tectonic plate identifier (default is "NOAM").
    date : int, optional
        ITRF reference year (default is 2014).
    grid_posting : float, optional
        Grid spacing in meters (default is 10000).

    Returns
    -------
    pmm_east_velocity : ndarray
        2D array of east-west velocity components (mm/year).
    pmm_north_velocity : ndarray
        2D array of north-south velocity components (mm/year).
    pmm_attributes : dict
        Dictionary containing grid metadata and model information.

    Raises
    ------
    ValueError
        If an unsupported ITRF year is specified.

    """
    # Validate ITRF date input
    _validate_itrf_year(date)

    # Generate regular grid for plate motion modeling
    logger.debug(f"Creating grid with {grid_posting/1000:.1f} km spacing")
    grid_longitudes, grid_latitudes = make_point_grid(frame_gdf, cell_size=grid_posting)

    # Extract grid dimensions for raster attributes
    grid_height, grid_width = grid_longitudes.shape

    # Create geospatial transform for the grid
    grid_transform = from_bounds(
        *frame_gdf.total_bounds, width=grid_width, height=grid_height
    )

    # Rigid-plate velocities (mm/yr) from the published ITRF pole (geepers)
    logger.info(f"Computing ITRF{date}-PMM plate velocities for plate '{plate}'")
    ve, vn, _ = plate_velocity_enu(
        grid_longitudes.ravel(), grid_latitudes.ravel(), plate, model=f"ITRF{date}-PMM"
    )
    pmm_east_velocity = ve.reshape((grid_height, grid_width))  # East-West component
    pmm_north_velocity = vn.reshape((grid_height, grid_width))  # North-South component

    # Compile grid and model metadata
    pmm_attributes = {
        "crs": frame_gdf.crs,
        "width": grid_width,
        "height": grid_height,
        "transform": grid_transform,
        "nodata": np.nan,
        "plate": plate,
        "model": f"ITRF{date}",
        "grid_posting_m": grid_posting,
        "units": "mm/year",
    }

    logger.info(f"Generated {grid_width}x{grid_height} velocity grid for {plate} plate")

    return pmm_east_velocity, pmm_north_velocity, pmm_attributes


def get_gia(
    frame_gdf: gpd.GeoDataFrame, model: str = "CARON"
) -> tuple[Any, dict[str, Any]]:
    """Load and process Glacial Isostatic Adjustment (GIA) data for a frame.

    Parameters
    ----------
    frame_gdf : GeoDataFrame
        GeoDataFrame representing the area of interest for clipping the GIA model.
    model : str, optional
        Name of the GIA model to use (default is 'CARON'). Options: 'CARON', 'ICE6'.

    Returns
    -------
    gia_raster : object
        Rasterized GIA data.
    gia_attr : dict
        Dictionary of attributes including the model name.

    Raises
    ------
    ValueError
        If an unsupported GIA model is specified.

    """
    # Normalize model name for comparison
    model_upper = model.upper()

    # Load the appropriate GIA model
    if model_upper == "CARON":
        gia_model = load_gia.load_caron_model(load_gia.CARON_GIA)

    elif model_upper == "ICE6":
        # Download and load ICE6G model data
        ice6_file_path = Path(load_gia.download_ice6g_data())
        gia_model = load_gia.load_ice6g_model(str(ice6_file_path))

        # Clean up temporary file
        with suppress(FileNotFoundError):
            ice6_file_path.unlink()

    else:
        msg = f"GIA model must be one of {SUPPORTED_GIA_MODELS}, got '{model}'"
        raise ValueError(msg)

    # Clip GIA model to the specified frame area of interest
    gia_frame = load_gia.clip_gia_df(gia_model, aoi=frame_gdf)

    # Rasterize the GIA data using vertical land motion rate
    # TODO: Add uncertainty
    gia_raster, gia_attributes = load_gia.rasterize_gdf(gia_frame, "vlm_rate")

    # Add model metadata to attributes
    gia_attributes["model"] = model

    return gia_raster, gia_attributes


def main(
    frame_id: int,
    output: str = "output",
    plate: Literal[
        "NOAM",
        "EURA",
        "PCFC",
        "NAZC",
        "SOAM",
        "AFRC",
        "ANTA",
        "ARAB",
        "INDI",
        "AUST",
        "CARI",
        "COCO",
        "JUAN",
        "PHIL",
    ] = "NOAM",
    itrf: Literal[2014, 2020] = 2014,
    gia: Literal["CARON", "ICE6"] = "CARON",
) -> None:
    """Generate DISP frame model with PMM and GIA data.

    Parameters
    ----------
    frame_id : int
        ID of the DISP frame to process.
    output : str, optional
        Output directory path (default is "output").
    plate : str, optional
        Tectonic plate identifier (default is "NOAM").
    itrf : {2014, 2020}, optional
        ITRF reference year (default is 2014).
    gia : {'CARON', 'ICE6'}, optional
        GIA model name (default is 'CARON').

    Raises
    ------
    FileNotFoundError
        If DISP frame database file is not found.
    ValueError
        If frame_id is not found in the database or invalid parameters are provided.

    """
    # Set output directory
    output_dir = Path(output)
    output_dir.mkdir(exist_ok=True)

    # Create output filename
    output_name = (
        output_dir / f"frame_{frame_id}_{plate}_ITRF{itrf}_GIA_{gia}_rates.tif"
    )

    # Get frame_df with opera_utils
    logger.info("Loading DISP frame database")
    disp_frame_db = opera_utils.get_frame_geodataframe()

    # Select requested frame
    selected_frame = disp_frame_db[disp_frame_db.index == frame_id]
    if selected_frame.empty:
        available_fids = disp_frame_db.index.unique()
        msg = f"Frame ID {frame_id} not found. Available IDs: {available_fids}"
        raise ValueError(msg)
    logger.info(f"Selected frame {frame_id}")

    # Get PMM for requested frame
    pmm_east_velocity, pmm_north_velocity, pmm_attributes = get_frame_pmm(
        selected_frame, plate=plate, date=itrf
    )

    # Get GIA for requested frame
    gia_model, gia_attributes = get_gia(selected_frame, model=gia)

    # Export to DISP frame
    output_frame = get_frame_utm_transform(selected_frame)
    resampling_mode = Resampling.bilinear

    with rasterio.open(
        output_name,
        "w",
        height=output_frame.height,
        width=output_frame.width,
        count=3,
        dtype="float32",
        crs=output_frame.crs,
        transform=output_frame.transform,
        nodata=np.nan,
        compress="lzw",
        tiled=True,
        blockxsize=512,
        blockysize=512,
        predictor=2,
        interleave="band",
    ) as dst:
        dst.set_band_description(1, "North-South (mm/yr)")
        dst.set_band_description(2, "East-West (mm/yr)")
        dst.set_band_description(3, "VLM (mm/yr)")

        # PMM reproject parameters
        pmm_reproject_params = {
            "src_transform": pmm_attributes["transform"],
            "src_crs": pmm_attributes["crs"],
            "dst_transform": output_frame.transform,
            "dst_crs": output_frame.crs,
            "resampling": resampling_mode,
            "dst_nodata": np.nan,
        }

        # Write PMM North-South rates
        logger.info(f"Writing plate NS rates: {output_name} [band: 1]")
        reproject(
            source=pmm_north_velocity,
            destination=rasterio.band(dst, 1),
            **pmm_reproject_params,
        )

        # Write PMM East-West rates
        logger.info(f"Writing plate EW rates: {output_name} [band: 2]")
        reproject(
            source=pmm_east_velocity,
            destination=rasterio.band(dst, 2),
            **pmm_reproject_params,
        )

        # GIA reproject parameters
        gia_reproject_params = {
            "src_transform": gia_attributes["transform"],
            "src_crs": gia_attributes["crs"],
            "dst_transform": output_frame.transform,
            "dst_crs": output_frame.crs,
            "resampling": resampling_mode,
            "dst_nodata": np.nan,
        }

        # Write GIA VLM rates
        logger.info(f"Writing GIA VLM rates: {output_name} [band: 3]")
        reproject(
            source=gia_model, destination=rasterio.band(dst, 3), **gia_reproject_params
        )

    logger.info(f"Successfully wrote output to: {output_name}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    tyro.cli(main)
