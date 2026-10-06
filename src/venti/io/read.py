# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Reading operations for geospatial raster data.

This module provides high-level classes and functions for reading
geospatial raster data with caching support.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class RasterMetadata:
    """Metadata for a raster dataset.

    Attributes
    ----------
    shape : tuple
        Raster shape (height, width) or (bands, height, width)
    transform : Any
        Affine transform
    crs : Any
        Coordinate reference system
    nodata : float, optional
        NoData value
    bounds : tuple, optional
        Geographic bounds (South, North, West, East)

    """

    shape: tuple[int, ...]
    transform: Any
    crs: Any
    nodata: float | None = None
    bounds: tuple[float, float, float, float] | None = None


@dataclass
class RasterData:
    """Container for raster data and metadata.

    Attributes
    ----------
    data : np.ndarray
        Raster data array
    metadata : RasterMetadata
        Raster metadata
    file_path : Path, optional
        Source file path

    """

    data: np.ndarray
    metadata: RasterMetadata
    file_path: Path | None = None

    @property
    def shape(self) -> tuple[int, ...]:
        """Get data shape."""
        return self.data.shape

    @property
    def transform(self) -> Any:
        """Get affine transform."""
        return self.metadata.transform

    @property
    def crs(self) -> Any:
        """Get CRS."""
        return self.metadata.crs

    @property
    def bounds(self) -> tuple[float, float, float, float] | None:
        """Get geographic bounds."""
        return self.metadata.bounds

    def write(
        self,
        output_path: str | Path,
        descriptions: list[str] | None = None,
    ) -> None:
        """Write raster data to file.

        Parameters
        ----------
        output_path : str or Path
            Output file path
        descriptions : list of str, optional
            Band descriptions

        """
        from .raster import write_geotiff

        write_geotiff(
            self.data,
            output_path,
            transform=self.metadata.transform,
            crs=self.metadata.crs,
            nodata=self.metadata.nodata,
            descriptions=descriptions,
        )
        logger.info(f"Written: {output_path}")


@dataclass
class NetCDFData:
    """Container for NetCDF data and metadata.

    Attributes
    ----------
    data : np.ndarray
        Data array
    coords : dict
        Coordinate arrays
    attrs : dict
        Attributes
    geo_info : dict
        Geographic information (CRS, transform, bounds)
    file_path : Path, optional
        Source file path

    """

    data: np.ndarray
    coords: dict[str, np.ndarray]
    attrs: dict[str, Any]
    geo_info: dict[str, Any]
    file_path: Path | None = None

    @property
    def shape(self) -> tuple[int, ...]:
        """Get data shape."""
        return self.data.shape

    @property
    def crs(self) -> Any:
        """Get CRS."""
        return self.geo_info.get("crs")

    @property
    def transform(self) -> Any:
        """Get affine transform."""
        return self.geo_info.get("transform")

    @property
    def bounds(self) -> tuple[float, float, float, float] | None:
        """Get geographic bounds."""
        return self.geo_info.get("bounds")


@dataclass
class RasterReader:
    """High-level interface for reading raster data.

    This class provides a clean interface for reading geospatial
    raster data in various formats with automatic caching.

    Attributes
    ----------
    cache : dict
        Cache for loaded rasters

    Examples
    --------
    ::

        reader = RasterReader()
        data = reader.read_geotiff('displacement.tif')
        print(f"Shape: {data.shape}, CRS: {data.crs}")

        # With caching
        data2 = reader.read_geotiff('displacement.tif')  # Uses cache
        print(f"Cache size: {reader.get_cache_size()}")

    """

    cache: dict[Path, RasterData] = field(default_factory=dict, init=False)

    def read_geotiff(
        self,
        file_path: str | Path,
        use_cache: bool = True,
    ) -> RasterData:
        """Read GeoTIFF file.

        Parameters
        ----------
        file_path : str or Path
            Path to GeoTIFF file
        use_cache : bool
            Use cached data if available

        Returns
        -------
        RasterData
            Raster data and metadata

        """
        from .raster import get_bounds, read_geotiff

        file_path = Path(file_path)

        if use_cache and file_path in self.cache:
            logger.debug(f"Using cached data: {file_path}")
            return self.cache[file_path]

        logger.info(f"Reading: {file_path}")
        data, geo_info = read_geotiff(file_path)

        # Get bounds
        try:
            bounds = get_bounds(file_path, as_latlon=False)
        except Exception:
            bounds = None

        metadata = RasterMetadata(
            shape=data.shape,
            transform=geo_info["transform"],
            crs=geo_info["crs"],
            nodata=geo_info.get("nodata"),
            bounds=bounds,
        )

        raster = RasterData(
            data=data,
            metadata=metadata,
            file_path=file_path,
        )

        if use_cache:
            self.cache[file_path] = raster

        return raster

    def read_netcdf(
        self,
        file_path: str | Path,
        variable: str = "displacement",
    ) -> NetCDFData:
        """Read NetCDF file.

        Parameters
        ----------
        file_path : str or Path
            Path to NetCDF file
        variable : str
            Variable name to read

        Returns
        -------
        NetCDFData
            NetCDF data and metadata

        """
        from .raster import read_netcdf

        file_path = Path(file_path)
        logger.info(f"Reading NetCDF: {file_path}")

        data, mask, geo_info = read_netcdf(file_path, variable=variable)

        return NetCDFData(
            data=data,
            coords={"mask": mask} if mask is not None else {},
            attrs={},  # Could be extended
            geo_info=geo_info,
            file_path=file_path,
        )

    def read_correction_layer(
        self,
        file_path: str | Path,
        layer: str,
    ) -> np.ndarray | None:
        """Read an OPERA DISP-S1 ``/corrections`` layer, if the product has it.

        Parameters
        ----------
        file_path : str or Path
            Path to the DISP-S1 NetCDF file.
        layer : str
            Correction layer name, e.g. ``'solid_earth_tide'``.

        Returns
        -------
        np.ndarray or None
            2D layer on the ``/displacement`` grid, or ``None`` if absent.

        """
        from .raster import read_netcdf_correction

        return read_netcdf_correction(file_path, layer)

    def get_bounds(
        self,
        file_path: str | Path,
        as_latlon: bool = False,
    ) -> tuple[float, float, float, float]:
        """Get bounds from a raster file.

        Parameters
        ----------
        file_path : str or Path
            Path to raster file
        as_latlon : bool
            Return bounds in lat/lon coordinates

        Returns
        -------
        tuple
            (South, North, West, East) bounds

        """
        from .raster import get_bounds

        return get_bounds(file_path, as_latlon=as_latlon)

    def clear_cache(self) -> None:
        """Clear the raster cache."""
        self.cache.clear()
        logger.debug("Cache cleared")

    def get_cache_size(self) -> int:
        """Get number of cached rasters.

        Returns
        -------
        int
            Number of cached rasters

        """
        return len(self.cache)
