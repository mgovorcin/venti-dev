"""I/O utilities for GeoTIFF and NetCDF files.

This module provides functions for reading and writing displacement data
in various formats with proper georeferencing.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

# Optional imports
try:
    import rasterio as rio
    from rasterio.transform import Affine, from_bounds

    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    Affine = None

try:
    import h5netcdf

    HAS_H5NETCDF = True
except ImportError:
    HAS_H5NETCDF = False


def read_geotiff(
    file_path: str | Path, band: int | None = None
) -> tuple[np.ndarray, dict[str, Any]]:
    """Read a GeoTIFF file and return data with georeferencing.

    Parameters
    ----------
    file_path : str or Path
        Path to GeoTIFF file
    band : int or None, optional
        Band number to read. If None, all bands are read and a 3D array
        (bands, rows, cols) is returned. For a single-band file with
        ``band=None``, the result is squeezed to 2D. Default is None.

    Returns
    -------
    data : np.ndarray
        2D array (single band) or 3D array (bands, rows, cols)
    geo_info : dict
        Dictionary with 'transform', 'crs', 'nodata'

    Raises
    ------
    ImportError
        If rasterio is not available

    Examples
    --------
    ::

        data, geo = read_geotiff('displacement.tif')
        print(f"CRS: {geo['crs']}")

    """
    if not HAS_RASTERIO:
        msg = "rasterio required. Install with: pip install rasterio"
        raise ImportError(msg)

    with rio.open(file_path) as src:
        if band is None:
            data = src.read()
            if data.shape[0] == 1:
                data = data.squeeze(axis=0)
        else:
            data = src.read(band)
        geo_info = {
            "transform": src.transform,
            "crs": src.crs,
            "nodata": src.nodata,
            "bounds": src.bounds,
            "shape": data.shape,
        }

    logger.info(f"Read GeoTIFF: {file_path}, shape={data.shape}")
    return data, geo_info


def write_geotiff(
    data: np.ndarray,
    output_path: str | Path,
    transform: Affine | None = None,
    crs: Any | None = None,
    nodata: float | None = None,
    reference_file: str | Path | None = None,
    dtype: str | None = None,
    compress: str = "lzw",
    descriptions: list[str] | None = None,
) -> None:
    """Write data to a GeoTIFF file.

    Parameters
    ----------
    data : np.ndarray
        2D or 3D array to write (3D: multiple bands)
    output_path : str or Path
        Output file path
    transform : Affine, optional
        Affine transform for georeferencing
    crs : Any, optional
        Coordinate reference system
    nodata : float, optional
        NoData value, by default None (uses np.nan for float types)
    reference_file : str or Path, optional
        Copy georeferencing from this file
    dtype : str, optional
        Output data type, by default None (use input dtype)
    compress : str, optional
        Compression method, by default 'lzw'
    descriptions : list of str, optional
        Band descriptions

    Raises
    ------
    ImportError
        If rasterio is not available
    ValueError
        If georeferencing information is missing

    Examples
    --------
    Write with explicit georeferencing::

        write_geotiff(
            corrected_disp,
            'corrected.tif',
            transform=transform,
            crs='EPSG:32610',
            nodata=np.nan
        )

    Write using reference file::

        write_geotiff(
            corrected_disp,
            'corrected.tif',
            reference_file='input.nc'
        )

    """
    if not HAS_RASTERIO:
        msg = "rasterio required. Install with: pip install rasterio"
        raise ImportError(msg)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get georeferencing from reference file if provided
    if reference_file is not None:
        ref_path = Path(reference_file)
        if ref_path.suffix.lower() == ".nc":
            _, _, geo_info = read_netcdf(reference_file)
            transform = geo_info.get("transform")
            crs = geo_info.get("crs")
            if nodata is None:
                nodata = geo_info.get("nodata")
        else:
            with rio.open(reference_file) as src:
                transform = src.transform
                crs = src.crs
                if nodata is None:
                    nodata = src.nodata

    # Validate georeferencing
    if transform is None or crs is None:
        msg = "Either (transform and crs) or reference_file must be provided"
        raise ValueError(msg)

    # Handle masked arrays
    if np.ma.isMaskedArray(data):
        if nodata is None:
            nodata = np.nan
        data = data.filled(nodata)

    # Set default nodata for float types
    if nodata is None and np.issubdtype(data.dtype, np.floating):
        nodata = np.nan

    # Determine output dtype
    if dtype is None:
        dtype = data.dtype

    # Handle dimensions
    if data.ndim == 2:
        height, width = data.shape
        count = 1
        data = data[np.newaxis, :, :]
    elif data.ndim == 3:
        count, height, width = data.shape
    else:
        msg = f"Data must be 2D or 3D, got {data.ndim}D"
        raise ValueError(msg)

    # Write file
    with rio.open(
        output_path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=count,
        dtype=dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
        compress=compress,
    ) as dst:
        dst.write(data)

        # Set band descriptions if provided
        if descriptions is not None:
            for i, desc in enumerate(descriptions, start=1):
                dst.set_band_description(i, desc)

    logger.info(f"Wrote GeoTIFF: {output_path}, shape={data.shape[1:]}, count={count}")


def read_netcdf(
    file_path: str | Path,
    variable: str = "displacement",
    mask_variable: str | None = "water_mask",
) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    """Read displacement data from NetCDF with georeferencing.

    Parameters
    ----------
    file_path : str or Path
        Path to NetCDF file
    variable : str, optional
        Variable name to read, by default 'displacement'
    mask_variable : str, optional
        Mask variable name, by default 'water_mask'
        Set to None to skip mask reading

    Returns
    -------
    data : np.ndarray
        2D array of data
    mask : np.ndarray or None
        2D mask array (if mask_variable provided)
    geo_info : dict
        Dictionary with 'transform', 'crs', 'nodata'

    Raises
    ------
    ValueError
        If required variables are missing

    Examples
    --------
    ::

        disp, mask, geo = read_netcdf('OPERA_displacement.nc')
        corrected = process(disp, mask)
        write_geotiff(corrected, 'output.tif', **geo)

    """
    ds = xr.open_dataset(file_path)

    # Read data variable
    if variable not in ds:
        msg = f"{variable.capitalize()} variable not found in {file_path}"
        raise ValueError(msg)

    data = ds[variable].values
    data = np.squeeze(data)

    if data.ndim != 2:
        msg = f"Expected 2D data, got {data.ndim}D"
        raise ValueError(msg)

    # Read mask if requested
    mask = None
    if mask_variable is not None:
        if mask_variable in ds:
            mask = ds[mask_variable].values
            mask = np.squeeze(mask)
            if mask.shape != data.shape:
                msg = f"Mask shape {mask.shape} doesn't match data shape {data.shape}"
                raise ValueError(msg)
        else:
            logger.warning(f"Mask variable '{mask_variable}' not found")

    # Extract georeferencing
    geo_info = _extract_netcdf_georef(ds, data.shape)

    # Get nodata
    data_var = ds[variable]
    if hasattr(data_var, "_FillValue"):
        geo_info["nodata"] = float(data_var._FillValue)
    elif hasattr(data_var, "missing_value"):
        geo_info["nodata"] = float(data_var.missing_value)
    else:
        geo_info["nodata"] = np.nan

    ds.close()

    logger.info(f"Read NetCDF: {file_path}, variable={variable}, shape={data.shape}")
    return data, mask, geo_info


def read_netcdf_correction(
    file_path: str | Path,
    layer: str,
    group: str = "corrections",
) -> np.ndarray | None:
    """Read a correction layer from an OPERA DISP-S1 NetCDF, if present.

    OPERA DISP-S1 products store correction layers (e.g.
    ``solid_earth_tide``, ``ionospheric_delay``) in a ``/corrections``
    group, on the same grid as ``/displacement``. Not every product
    version carries every layer, so a missing group or layer is an
    expected case and returns ``None`` instead of raising.

    Parameters
    ----------
    file_path : str or Path
        Path to the DISP-S1 NetCDF file.
    layer : str
        Variable name inside `group`, e.g. ``'solid_earth_tide'``.
    group : str, optional
        NetCDF group holding the layer, by default ``'corrections'``.

    Returns
    -------
    np.ndarray or None
        2D correction layer, or ``None`` if `group` or `layer` is absent.

    Raises
    ------
    ValueError
        If the layer is present but not 2D.

    Examples
    --------
    ::

        set_corr = read_netcdf_correction(disp_file, "solid_earth_tide")
        if set_corr is not None:
            disp = disp - set_corr

    """
    # Only xarray opens the file: its reads share a lock around HDF5, which the
    # pip netCDF4 wheels need when `CalibrationWorkflow.run` uses threads.
    try:
        ds = xr.open_dataset(file_path, group=group)
    except OSError as exc:
        if "group not found" not in str(exc):
            raise
        return None
    with ds:
        if layer not in ds.variables:
            return None
        data = np.squeeze(ds[layer].values)

    if data.ndim != 2:
        msg = f"Expected 2D {group}/{layer} in {file_path}, got {data.ndim}D"
        raise ValueError(msg)
    return data


def update_netcdf_variable(
    file_path: str | Path,
    variable_name: str,
    data: np.ndarray,
    variable_attrs: dict[str, Any] | None = None,
    create_if_missing: bool = True,
) -> None:
    """Update or create a variable in an existing NetCDF file.

    Parameters
    ----------
    file_path : str or Path
        Path to NetCDF file
    variable_name : str
        Name of variable to update/create
    data : np.ndarray
        Data to write
    variable_attrs : dict, optional
        Variable attributes to set
    create_if_missing : bool, optional
        Create variable if it doesn't exist, by default True

    Raises
    ------
    ImportError
        If h5netcdf is not available

    Examples
    --------
    ::

        update_netcdf_variable(
            'displacement.nc',
            'displacement_corrected',
            corrected_data,
            variable_attrs={'units': 'm', 'description': 'Corrected displacement'}
        )

    """
    if not HAS_H5NETCDF:
        msg = "h5netcdf required. Install with: pip install h5netcdf"
        raise ImportError(msg)

    def _handle_missing_variable() -> None:
        """Handle missing variable error."""
        if create_if_missing:
            logger.warning(
                f"Variable '{variable_name}' not found. Creating new variables in"
                " NetCDF requires careful dimension handling."
            )
            msg = (
                "Creating new NetCDF variables not yet implemented. "
                "Use write_geotiff or xarray.to_netcdf instead."
            )
            raise NotImplementedError(msg)
        msg = f"Variable '{variable_name}' not found in {file_path}"
        raise ValueError(msg)

    try:
        with h5netcdf.File(file_path, "r+") as f:
            if variable_name in f.variables:
                # Update existing variable
                var = f.variables[variable_name]
                var[...] = data
                logger.info(f"Updated variable '{variable_name}' in {file_path}")
            else:
                _handle_missing_variable()

            # Update attributes if provided
            if variable_attrs is not None and variable_name in f.variables:
                var = f.variables[variable_name]
                var.attrs.update(variable_attrs)

    except Exception:
        logger.exception("Error updating NetCDF")
        raise


def get_bounds(
    file_path: str | Path, as_latlon: bool = False
) -> tuple[float, float, float, float]:
    """Get bounding box from raster file.

    Parameters
    ----------
    file_path : str or Path
        Path to raster file (GeoTIFF or NetCDF)
    as_latlon : bool, optional
        Return bounds in lat/lon (EPSG:4326), by default False

    Returns
    -------
    tuple
        Bounds as (South, North, West, East)

    Examples
    --------
    ::

        bounds = get_bounds('displacement.nc')
        S, N, W, E = bounds

    """
    file_path = Path(file_path)

    if file_path.suffix.lower() == ".nc":
        # NetCDF file - use xarray
        with xr.open_dataset(file_path) as ds:
            if "x" in ds.coords and "y" in ds.coords:
                x = ds.coords["x"].values
                y = ds.coords["y"].values
                W, E = float(x.min()), float(x.max())
                S, N = float(y.min()), float(y.max())
            else:
                msg = "No x/y coordinates found in NetCDF"
                raise ValueError(msg)

            # Transform to lat/lon if requested
            if as_latlon and HAS_RASTERIO:
                # Try to get CRS
                crs = _extract_netcdf_georef(ds, (len(y), len(x)))["crs"]
                if crs and crs != "EPSG:4326":
                    bounds = rio.warp.transform_bounds(
                        crs, "EPSG:4326", W, S, E, N, densify_pts=21
                    )
                    W, S, E, N = bounds[0], bounds[1], bounds[2], bounds[3]

    else:
        # Assume GeoTIFF
        if not HAS_RASTERIO:
            msg = "rasterio required for GeoTIFF files"
            raise ImportError(msg)

        with rio.open(file_path) as src:
            bounds = src.bounds
            W, S, E, N = bounds.left, bounds.bottom, bounds.right, bounds.top

            if as_latlon and src.crs != "EPSG:4326":
                bounds = rio.warp.transform_bounds(
                    src.crs, "EPSG:4326", W, S, E, N, densify_pts=21
                )
                W, S, E, N = bounds[0], bounds[1], bounds[2], bounds[3]

    return S, N, W, E


def _extract_netcdf_georef(
    ds: xr.Dataset, _data_shape: tuple[int, int]
) -> dict[str, Any]:
    """Extract georeferencing from NetCDF dataset."""
    geo_info = {}

    # Get CRS
    if "spatial_ref" in ds:
        sr = ds["spatial_ref"]
        for crs_attr in ["crs_wkt", "spatial_ref", "wkt"]:
            if crs_attr in sr.attrs:
                geo_info["crs"] = sr.attrs[crs_attr]
                break
    elif "crs" in ds.attrs:
        geo_info["crs"] = ds.attrs["crs"]
    else:
        geo_info["crs"] = "EPSG:4326"
        logger.warning("No CRS found, assuming EPSG:4326")

    # Get transform
    if not HAS_RASTERIO:
        logger.warning("rasterio not available, skipping transform extraction")
        return geo_info

    if "spatial_ref" in ds and "GeoTransform" in ds["spatial_ref"].attrs:
        # Parse GeoTransform
        gt_str = ds["spatial_ref"].attrs["GeoTransform"]
        gt = [float(x) for x in gt_str.split()]
        transform = Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3])
        geo_info["transform"] = transform

    elif "x" in ds.coords and "y" in ds.coords:
        # Build from coordinates
        x = ds.coords["x"].values
        y = ds.coords["y"].values

        x_res = abs(x[1] - x[0]) if len(x) > 1 else 1.0
        y_res = abs(y[1] - y[0]) if len(y) > 1 else 1.0

        x_min = float(x.min()) - x_res / 2
        x_max = float(x.max()) + x_res / 2
        y_min = float(y.min()) - y_res / 2
        y_max = float(y.max()) + y_res / 2

        transform = from_bounds(x_min, y_min, x_max, y_max, len(x), len(y))
        geo_info["transform"] = transform

    return geo_info
