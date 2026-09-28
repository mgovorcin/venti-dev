"""Spatial tools for InSAR calibration.

Windowed surface fitting (`SpatialProcessor`), interpolation of scattered
points onto a raster grid, gap filling (`fill_gaps`: inverse-distance, used
before the fit; `fill_masked_region`: nearest-neighbour, used for event
masks), residual-region detection and resampling.
"""

from __future__ import annotations

from .gap_filling import fill_gaps
from .interpolation import (
    detect_coherent_residual_regions,
    fill_masked_region,
    grid_coordinates,
    interpolate_griddata,
    interpolate_rbf,
)
from .processor import SpatialProcessor
from .resample import downsample_array, upsample_array

__all__ = [
    "SpatialProcessor",
    "detect_coherent_residual_regions",
    "downsample_array",
    "fill_gaps",
    "fill_masked_region",
    "grid_coordinates",
    "interpolate_griddata",
    "interpolate_rbf",
    "upsample_array",
]
