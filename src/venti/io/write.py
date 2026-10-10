# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Writing operations for geospatial raster data.

This module provides high-level classes and functions for writing
geospatial raster data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class RasterWriter:
    """High-level interface for writing raster data.

    This class provides a clean interface for writing geospatial
    raster data in various formats.

    Examples
    --------
    ::

        writer = RasterWriter()
        writer.write_geotiff(
            data,
            'output.tif',
            transform=transform,
            crs='EPSG:4326'
        )

    """

    def write_geotiff(
        self,
        data: np.ndarray,
        output_path: str | Path,
        reference_file: str | Path | None = None,
        transform: Any | None = None,
        crs: Any | None = None,
        nodata: float | None = None,
        descriptions: list[str] | None = None,
    ) -> None:
        """Write GeoTIFF file.

        Parameters
        ----------
        data : np.ndarray
            Data to write
        output_path : str or Path
            Output file path
        reference_file : str or Path, optional
            Reference file for georeferencing
        transform : Any, optional
            Affine transform
        crs : Any, optional
            Coordinate reference system
        nodata : float, optional
            NoData value
        descriptions : list of str, optional
            Band descriptions

        """
        from .raster import write_geotiff

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        write_geotiff(
            data,
            output_path,
            reference_file=reference_file,
            transform=transform,
            crs=crs,
            nodata=nodata,
            descriptions=descriptions,
        )
        logger.info(f"Written: {output_path}")
