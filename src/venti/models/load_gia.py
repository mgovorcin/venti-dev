# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""GIA (Glacial Isostatic Adjustment) data processing utilities.

This module provides functions to download, load, and process GIA model data
from various sources including ICE6G-D and Caron et al. 2018 models.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests  # type: ignore[import-untyped]
import xarray as xr
from rasterio.transform import from_bounds
from scipy.interpolate import griddata

# Constants
ICE6D_URL = (
    "https://www.atmosp.physics.utoronto.ca/~peltier/datasets/"
    "Ice6G_D_VM5a_O512/drad.12mgrid_512.nc"
)
# Caron et al 2018, https://doi.org/10.1002/2017GL076644
CARON_GIA = Path(__file__).parent / "data/GIA_maps_Caron_et_al_2018"
DEFAULT_BUFFER_DISTANCE = 300e3  # 300 km in meters


class GIALoadError(Exception):
    """Raised when GIA model download or load fails."""


def download_ice6g_data(url=ICE6D_URL, local_filename="ice6g_data.nc"):
    """Download ICE6G-D Glacial Isostatic Adjustment (GIA) model data.

    The ICE6G-D model provides estimates of glacial isostatic adjustment
    (GIA) velocities and uplift rates on a global grid. This function
    downloads the NetCDF file from a remote repository if it is not already
    present in the specified directory.

    Parameters
    ----------
    url : str, optional
        The URL where the ICE6G-D NetCDF file can be downloaded.
        Defaults to the official ICE6G-D data source.
    local_filename : str, optional
        The name of the NetCDF file to save locally. Defaults to
        ``"ice6g_data.nc"``.

    Returns
    -------
    Path
        Path to the downloaded NetCDF file.

    Raises
    ------
    GIADataDownloadError
        If the download fails or the remote file is not reachable.

    """
    # Check if file already exists
    local_path = Path(local_filename)

    # Check if file already exists
    if local_path.exists():
        print(f"File {local_filename} already exists, loading...")
        return local_path

    # Download the file
    print(f"Downloading from {url}...")
    response = requests.get(url)

    if response.status_code == 200:
        print(f"Download successful! ({len(response.content)} bytes)")

        # Save to file
        with open(local_filename, "wb") as f:
            f.write(response.content)

        return Path(local_filename).resolve()
    else:
        msg = f"Download failed: {response.status_code}"
        raise GIALoadError(msg)


def load_ice6g_model(file_path: str | Path) -> gpd.GeoDataFrame:
    """Load ICE6G-D GIA model data and return a GeoDataFrame.

    The ICE6G-D NetCDF file is parsed to extract latitude, longitude,
    and vertical velocity fields, which are converted to a GeoDataFrame
    in WGS84 coordinates.

    Reference: Peltier et al 2018, doi:10.1002/2016JB013844


    Parameters
    ----------
    file_path : str
        Path to the ICE6G-D NetCDF file.

    Returns
    -------
    geopandas.GeoDataFrame
        A GeoDataFrame containing point geometries and vertical
        velocity estimates.

    Raises
    ------
    FileNotFoundError
        If the input NetCDF file does not exist.

    """
    if file_path is None:
        msg = "ice6g_path must be provided when loading ICE6G-D model"
        raise ValueError(msg)

    # Read NetCDF data
    da_ice6g = xr.open_dataset(file_path)
    ice6_df = da_ice6g.to_dataframe().reset_index()

    # Shift longitude from 0-360 to -180-180
    ice6_df.loc[ice6_df.Lon > 180, "Lon"] = ice6_df.loc[ice6_df.Lon > 180, "Lon"] - 360

    # Create GeoDataFrame
    ice6_gdf = gpd.GeoDataFrame(
        ice6_df, geometry=gpd.points_from_xy(ice6_df.Lon, ice6_df.Lat), crs="EPSG:4326"
    )

    # Rename columns for consistency
    ice6_gdf = ice6_gdf.rename(
        columns={"Drad_250": "vlm_rate", "Lon": "lon", "Lat": "lat"}
    )
    return ice6_gdf


def load_caron_model(file_path: str | Path) -> gpd.GeoDataFrame:
    """Load and process Caron et al. 2018 GIA model."""
    if file_path is None:
        msg = "caron_path must be provided when loading Caron model"
        raise ValueError(msg)

    # Define column names
    names = [
        "lat",
        "lon",
        "vlm_rate",
        "vlm_std",
        "geoid_rate",
        "geoid_std",
        "gravity_rate",
        "gravity_std",
    ]

    # Read the data
    gia_df = pd.read_csv(
        file_path, skiprows=6, delimiter=r"\s+", names=names, dtype=float
    )

    # Shift longitude from 0-360 to -180-180
    gia_df.loc[gia_df.lon > 180, "lon"] = gia_df.loc[gia_df.lon > 180, "lon"] - 360

    # Adjust latitude coordinates
    gia_df.lat -= 90
    gia_df.lat = np.flipud(gia_df.lat)

    # Create GeoDataFrame
    caron_gdf = gpd.GeoDataFrame(
        gia_df, geometry=gpd.points_from_xy(gia_df.lon, gia_df.lat), crs="EPSG:4326"
    )
    return caron_gdf


def clip_gia_df(
    gdf: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame,
    buffer_distance: float = DEFAULT_BUFFER_DISTANCE,
) -> gpd.GeoDataFrame:
    """Clip GIA data to an area of interest (AOI) with an optional buffer.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Input GeoDataFrame containing GIA data points.
    aoi : geopandas.GeoDataFrame
        Polygon or multipolygon defining the area of interest.
    buffer_distance : float, optional
        Buffer distance around the AOI in kilometers. Default is 50 km.

    Returns
    -------
    geopandas.GeoDataFrame
        Clipped subset of the input GIA dataset.

    """
    # Get aoi in utm
    utm_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(utm_crs)

    # Add buffer to aoi and reproject to gdf.crs
    aoi_expanded = aoi_utm.buffer(buffer_distance).to_crs(gdf.crs)

    # Clip
    gdf_clipped = gdf.clip(aoi_expanded)

    return gdf_clipped


def rasterize_gdf(gdf: gpd.GeoDataFrame, value_column: str) -> tuple[np.array, dict]:
    """Rasterize GIA point data to a regular grid.

    The function interpolates vertical velocity values from a GeoDataFrame
    of GIA points onto a regular latitude longitude grid using inverse
    distance weighting.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Input GIA data with geometry and value fields.
    value_column : str
        Column name in the GeoDataFrame containing the values to interpolate.

    Returns
    -------
    tuple of (numpy.ndarray, dict)
        A tuple containing the rasterized data array and associated metadata.

    """
    # Get bounds from the dataframe
    minx, miny, maxx, maxy = gdf.total_bounds

    # Define grid resolution (adjust as needed)
    unique_lons = np.sort(gdf.lon.unique())
    unique_lats = np.sort(gdf.lat.unique())

    width = len(unique_lons)
    height = len(unique_lats)

    # Create the transform
    src_transform = from_bounds(minx, miny, maxx, maxy, width, height)

    # Create coordinate grids directly
    grid_lons, grid_lats = np.meshgrid(unique_lons, unique_lats)

    # Extract point coordinates and values from the dataframe
    coords = np.column_stack((gdf.geometry.x, gdf.geometry.y))
    values = gdf[value_column].values

    # Perform gridded interpolation
    grid_vals = griddata(
        points=coords,  # (N,2) array of point coordinates
        values=values,  # (N,) values
        xi=(grid_lons, grid_lats),  # 2D grid coordinates
        method="linear",  # options: 'linear', 'nearest', 'cubic'
        fill_value=np.nan,  # Fill value for points outside convex hull
    )

    grid_attr = {
        "crs": gdf.crs,
        "width": width,
        "height": height,
        "transform": src_transform,
        "nodata": np.nan,
    }
    return grid_vals, grid_attr
