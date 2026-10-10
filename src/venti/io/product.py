# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Product classes for creating Venti NetCDF output products.

This module provides classes for creating standardized NetCDF products:
- CalProduct: Calibration products with GNSS-corrected displacement
- VlmProduct: Decomposition products with East-North-Up components
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

# Type checking imports
if TYPE_CHECKING:
    from rasterio.transform import Affine


@dataclass
class CalProduct:
    """Calibration product for GNSS-corrected InSAR displacement.

    Creates a NetCDF product containing:
    - Calibrated displacement (corrected with GNSS)
    - GNSS east, north, up components used for calibration
    - Standard deviations for all components
    - Metadata and georeferencing information

    Attributes
    ----------
    output_path : Path
        Path to output NetCDF file
    calibrated_displacement : np.ndarray
        2D array of calibrated displacement values
    gnss_east : np.ndarray
        2D array of GNSS east component
    gnss_north : np.ndarray
        2D array of GNSS north component
    gnss_up : np.ndarray
        2D array of GNSS up component
    displacement_std : np.ndarray, optional
        Standard deviation of displacement
    gnss_east_std : np.ndarray, optional
        Standard deviation of GNSS east component
    gnss_north_std : np.ndarray, optional
        Standard deviation of GNSS north component
    gnss_up_std : np.ndarray, optional
        Standard deviation of GNSS up component
    transform : Affine, optional
        Affine transform for georeferencing
    crs : Any, optional
        Coordinate reference system
    x_coords : np.ndarray, optional
        X coordinate array
    y_coords : np.ndarray, optional
        Y coordinate array
    reference_datetime : datetime, optional
        Reference date for displacement
    secondary_datetime : datetime, optional
        Secondary date for displacement
    metadata : dict, optional
        Additional metadata to include in product

    Examples
    --------
    ::

        product = CalProduct(
            output_path='calibrated_displacement.nc',
            calibrated_displacement=disp_data,
            gnss_east=gnss_e,
            gnss_north=gnss_n,
            gnss_up=gnss_u,
            displacement_std=disp_std,
            transform=transform,
            crs='EPSG:32610'
        )
        product.write()

    """

    output_path: Path
    calibrated_displacement: np.ndarray
    gnss_east: np.ndarray
    gnss_north: np.ndarray
    gnss_up: np.ndarray
    displacement_std: np.ndarray | None = None
    gnss_east_std: np.ndarray | None = None
    gnss_north_std: np.ndarray | None = None
    gnss_up_std: np.ndarray | None = None
    transform: Affine | None = None
    crs: Any | None = None
    x_coords: np.ndarray | None = None
    y_coords: np.ndarray | None = None
    reference_datetime: datetime | None = None
    secondary_datetime: datetime | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate inputs and convert path."""
        self.output_path = Path(self.output_path)

        # Validate array shapes
        shape = self.calibrated_displacement.shape
        arrays_to_check = [
            ("gnss_east", self.gnss_east),
            ("gnss_north", self.gnss_north),
            ("gnss_up", self.gnss_up),
        ]

        for name, array in arrays_to_check:
            if array.shape != shape:
                msg = (
                    f"{name} shape {array.shape} does not match "
                    f"displacement shape {shape}"
                )
                raise ValueError(msg)

        # Validate standard deviation arrays if provided
        std_arrays = [
            ("displacement_std", self.displacement_std),
            ("gnss_east_std", self.gnss_east_std),
            ("gnss_north_std", self.gnss_north_std),
            ("gnss_up_std", self.gnss_up_std),
        ]

        for name, array in std_arrays:
            if array is not None and array.shape != shape:
                msg = f"{name} shape {array.shape} does not match shape {shape}"
                raise ValueError(msg)

        # Generate coordinates if not provided
        if self.x_coords is None or self.y_coords is None:
            self._generate_coordinates()

    def _generate_coordinates(self):
        """Generate coordinate arrays from transform."""
        if self.transform is None:
            # Create simple pixel coordinates
            rows, cols = self.calibrated_displacement.shape
            self.y_coords = np.arange(rows, dtype=np.float64)
            self.x_coords = np.arange(cols, dtype=np.float64)
            logger.warning("No transform provided, using pixel coordinates")
        else:
            # Generate real-world coordinates from affine transform
            rows, cols = self.calibrated_displacement.shape
            x_pixel = np.arange(cols)
            y_pixel = np.arange(rows)

            # Transform to real-world coordinates
            self.x_coords = self.transform[0] + x_pixel * self.transform[1]
            self.y_coords = self.transform[3] + y_pixel * self.transform[5]

    def to_dataset(self) -> xr.Dataset:
        """Create xarray Dataset from product data.

        Returns
        -------
        xr.Dataset
            Dataset containing all product layers with metadata

        """
        # Create coordinate arrays
        coords = {
            "y": self.y_coords,
            "x": self.x_coords,
        }

        # Create data variables
        data_vars = {
            "calibrated_displacement": (
                ["y", "x"],
                self.calibrated_displacement,
                {
                    "units": "meters",
                    "long_name": "Calibrated Line-of-Sight Displacement",
                    "description": (
                        "InSAR displacement corrected using GNSS reference data"
                    ),
                },
            ),
            "gnss_east": (
                ["y", "x"],
                self.gnss_east,
                {
                    "units": "meters",
                    "long_name": "GNSS East Component",
                    "description": "East component of GNSS displacement field",
                },
            ),
            "gnss_north": (
                ["y", "x"],
                self.gnss_north,
                {
                    "units": "meters",
                    "long_name": "GNSS North Component",
                    "description": "North component of GNSS displacement field",
                },
            ),
            "gnss_up": (
                ["y", "x"],
                self.gnss_up,
                {
                    "units": "meters",
                    "long_name": "GNSS Up Component",
                    "description": "Vertical (up) component of GNSS displacement field",
                },
            ),
        }

        # Add standard deviation layers if available
        if self.displacement_std is not None:
            data_vars["displacement_std"] = (
                ["y", "x"],
                self.displacement_std,
                {
                    "units": "meters",
                    "long_name": "Displacement Standard Deviation",
                    "description": "Uncertainty in calibrated displacement",
                },
            )

        if self.gnss_east_std is not None:
            data_vars["gnss_east_std"] = (
                ["y", "x"],
                self.gnss_east_std,
                {
                    "units": "meters",
                    "long_name": "GNSS East Standard Deviation",
                    "description": "Uncertainty in GNSS east component",
                },
            )

        if self.gnss_north_std is not None:
            data_vars["gnss_north_std"] = (
                ["y", "x"],
                self.gnss_north_std,
                {
                    "units": "meters",
                    "long_name": "GNSS North Standard Deviation",
                    "description": "Uncertainty in GNSS north component",
                },
            )

        if self.gnss_up_std is not None:
            data_vars["gnss_up_std"] = (
                ["y", "x"],
                self.gnss_up_std,
                {
                    "units": "meters",
                    "long_name": "GNSS Up Standard Deviation",
                    "description": "Uncertainty in GNSS up component",
                },
            )

        # Create dataset
        ds = xr.Dataset(data_vars=data_vars, coords=coords)

        # Add global attributes
        attrs: dict[str, Any] = {
            "title": "Venti Calibration Product",
            "product_type": "CAL",
            "institution": "NASA JPL",
            "source": "Venti - GNSS-InSAR Calibration Toolkit",
            "history": f"Created {datetime.now().isoformat()}",
        }

        # Add CRS if available
        if self.crs is not None:
            attrs["crs"] = str(self.crs)

        # Add date information
        if self.reference_datetime is not None:
            attrs["reference_datetime"] = self.reference_datetime.isoformat()
        if self.secondary_datetime is not None:
            attrs["secondary_datetime"] = self.secondary_datetime.isoformat()

        # Add affine transform as attributes
        if self.transform is not None:
            attrs["geotransform"] = list(self.transform)[:6]

        # Add custom metadata
        attrs.update(self.metadata)

        ds.attrs = attrs

        return ds

    def write(
        self,
        compression: str = "zlib",
        complevel: int = 4,
        engine: str = "h5netcdf",
    ) -> None:
        """Write product to NetCDF file.

        Parameters
        ----------
        compression : str, optional
            Compression algorithm, by default 'zlib'
        complevel : int, optional
            Compression level (1-9), by default 4
        engine : str, optional
            NetCDF engine to use, by default 'h5netcdf'

        """
        # Create output directory if needed
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

        # Create dataset
        ds = self.to_dataset()

        # Encoding for compression
        encoding = {}
        for var in ds.data_vars:
            encoding[var] = {
                "zlib": compression == "zlib",
                "complevel": complevel,
                "dtype": "float32",
            }

        # Write to file
        logger.info(f"Writing calibration product to {self.output_path}")
        ds.to_netcdf(
            self.output_path,
            mode="w",
            format="NETCDF4",
            engine=engine,
            encoding=encoding,
        )
        logger.info(
            f"Created calibration product: {self.output_path} "
            f"({self.output_path.stat().st_size / 1e6:.2f} MB)"
        )


@dataclass
class VlmProduct:
    """Vertical Land Motion (VLM) product with decomposed displacement.

    Creates a NetCDF product containing:
    - East, North, Up displacement components from InSAR decomposition
    - Standard deviations for all components
    - Metadata and georeferencing information

    Attributes
    ----------
    output_path : Path
        Path to output NetCDF file
    east_displacement : np.ndarray
        2D array of east displacement component
    north_displacement : np.ndarray
        2D array of north displacement component
    up_displacement : np.ndarray
        2D array of up (vertical) displacement component
    east_std : np.ndarray, optional
        Standard deviation of east component
    north_std : np.ndarray, optional
        Standard deviation of north component
    up_std : np.ndarray, optional
        Standard deviation of up component
    transform : Affine, optional
        Affine transform for georeferencing
    crs : Any, optional
        Coordinate reference system
    x_coords : np.ndarray, optional
        X coordinate array
    y_coords : np.ndarray, optional
        Y coordinate array
    reference_datetime : datetime, optional
        Reference date for displacement
    secondary_datetime : datetime, optional
        Secondary date for displacement
    metadata : dict, optional
        Additional metadata to include in product

    Examples
    --------
    ::

        product = VlmProduct(
            output_path='decomposed_displacement.nc',
            east_displacement=east_data,
            north_displacement=north_data,
            up_displacement=up_data,
            east_std=east_std,
            north_std=north_std,
            up_std=up_std,
            transform=transform,
            crs='EPSG:32610'
        )
        product.write()

    """

    output_path: Path
    east_displacement: np.ndarray
    north_displacement: np.ndarray
    up_displacement: np.ndarray
    east_std: np.ndarray | None = None
    north_std: np.ndarray | None = None
    up_std: np.ndarray | None = None
    transform: Affine | None = None
    crs: Any | None = None
    x_coords: np.ndarray | None = None
    y_coords: np.ndarray | None = None
    reference_datetime: datetime | None = None
    secondary_datetime: datetime | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate inputs and convert path."""
        self.output_path = Path(self.output_path)

        # Validate array shapes
        shape = self.east_displacement.shape
        arrays_to_check = [
            ("north_displacement", self.north_displacement),
            ("up_displacement", self.up_displacement),
        ]

        for name, array in arrays_to_check:
            if array.shape != shape:
                msg = (
                    f"{name} shape {array.shape} does not match "
                    f"east_displacement shape {shape}"
                )
                raise ValueError(msg)

        # Validate standard deviation arrays if provided
        std_arrays = [
            ("east_std", self.east_std),
            ("north_std", self.north_std),
            ("up_std", self.up_std),
        ]

        for name, array in std_arrays:
            if array is not None and array.shape != shape:
                msg = f"{name} shape {array.shape} does not match shape {shape}"
                raise ValueError(msg)

        # Generate coordinates if not provided
        if self.x_coords is None or self.y_coords is None:
            self._generate_coordinates()

    def _generate_coordinates(self):
        """Generate coordinate arrays from transform."""
        if self.transform is None:
            # Create simple pixel coordinates
            rows, cols = self.east_displacement.shape
            self.y_coords = np.arange(rows, dtype=np.float64)
            self.x_coords = np.arange(cols, dtype=np.float64)
            logger.warning("No transform provided, using pixel coordinates")
        else:
            # Generate real-world coordinates from affine transform
            rows, cols = self.east_displacement.shape
            x_pixel = np.arange(cols)
            y_pixel = np.arange(rows)

            # Transform to real-world coordinates
            self.x_coords = self.transform[0] + x_pixel * self.transform[1]
            self.y_coords = self.transform[3] + y_pixel * self.transform[5]

    def to_dataset(self) -> xr.Dataset:
        """Create xarray Dataset from product data.

        Returns
        -------
        xr.Dataset
            Dataset containing all product layers with metadata

        """
        # Create coordinate arrays
        coords = {
            "y": self.y_coords,
            "x": self.x_coords,
        }

        # Create data variables
        data_vars = {
            "east": (
                ["y", "x"],
                self.east_displacement,
                {
                    "units": "meters",
                    "long_name": "East Displacement Component",
                    "description": (
                        "East component of 3D displacement from InSAR decomposition"
                    ),
                },
            ),
            "north": (
                ["y", "x"],
                self.north_displacement,
                {
                    "units": "meters",
                    "long_name": "North Displacement Component",
                    "description": (
                        "North component of 3D displacement from InSAR decomposition"
                    ),
                },
            ),
            "up": (
                ["y", "x"],
                self.up_displacement,
                {
                    "units": "meters",
                    "long_name": "Vertical Displacement Component",
                    "description": (
                        "Vertical (up) component of 3D displacement "
                        "from InSAR decomposition"
                    ),
                },
            ),
        }

        # Add standard deviation layers if available
        if self.east_std is not None:
            data_vars["east_std"] = (
                ["y", "x"],
                self.east_std,
                {
                    "units": "meters",
                    "long_name": "East Component Standard Deviation",
                    "description": "Uncertainty in east displacement component",
                },
            )

        if self.north_std is not None:
            data_vars["north_std"] = (
                ["y", "x"],
                self.north_std,
                {
                    "units": "meters",
                    "long_name": "North Component Standard Deviation",
                    "description": "Uncertainty in north displacement component",
                },
            )

        if self.up_std is not None:
            data_vars["up_std"] = (
                ["y", "x"],
                self.up_std,
                {
                    "units": "meters",
                    "long_name": "Vertical Component Standard Deviation",
                    "description": "Uncertainty in vertical displacement component",
                },
            )

        # Create dataset
        ds = xr.Dataset(data_vars=data_vars, coords=coords)

        # Add global attributes
        attrs: dict[str, Any] = {
            "title": "Venti VLM Product",
            "product_type": "VLM",
            "institution": "NASA JPL",
            "source": "Venti - InSAR Displacement Decomposition",
            "history": f"Created {datetime.now().isoformat()}",
        }

        # Add CRS if available
        if self.crs is not None:
            attrs["crs"] = str(self.crs)

        # Add date information
        if self.reference_datetime is not None:
            attrs["reference_datetime"] = self.reference_datetime.isoformat()
        if self.secondary_datetime is not None:
            attrs["secondary_datetime"] = self.secondary_datetime.isoformat()

        # Add affine transform as attributes
        if self.transform is not None:
            attrs["geotransform"] = list(self.transform)[:6]

        # Add custom metadata
        attrs.update(self.metadata)

        ds.attrs = attrs

        return ds

    def write(
        self,
        compression: str = "zlib",
        complevel: int = 4,
        engine: str = "h5netcdf",
    ) -> None:
        """Write product to NetCDF file.

        Parameters
        ----------
        compression : str, optional
            Compression algorithm, by default 'zlib'
        complevel : int, optional
            Compression level (1-9), by default 4
        engine : str, optional
            NetCDF engine to use, by default 'h5netcdf'

        """
        # Create output directory if needed
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

        # Create dataset
        ds = self.to_dataset()

        # Encoding for compression
        encoding = {}
        for var in ds.data_vars:
            encoding[var] = {
                "zlib": compression == "zlib",
                "complevel": complevel,
                "dtype": "float32",
            }

        # Write to file
        logger.info(f"Writing VLM product to {self.output_path}")
        ds.to_netcdf(
            self.output_path,
            mode="w",
            format="NETCDF4",
            engine=engine,
            encoding=encoding,
        )
        logger.info(
            f"Created VLM product: {self.output_path} "
            f"({self.output_path.stat().st_size / 1e6:.2f} MB)"
        )
