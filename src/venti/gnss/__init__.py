"""GNSS reference fields for InSAR calibration.

Download UNR gridded GNSS products (`GNSSReference`, or the functions in
`venti.gnss.unr`), read station rates or epoch displacements, and project
them into an InSAR line-of-sight field on any raster grid.
"""

from __future__ import annotations

from .los import project_to_los, project_uncertainty_to_los
from .reference import GNSSReference, compute_gnss_los, compute_gnss_los_std
from .unr import (
    calculate_station_velocity,
    download_grid_lookup,
    download_station,
    find_stations_in_bounds,
    read_epoch_displacements,
    read_station_rate,
)

__all__ = [
    "GNSSReference",
    "calculate_station_velocity",
    "compute_gnss_los",
    "compute_gnss_los_std",
    "download_grid_lookup",
    "download_station",
    "find_stations_in_bounds",
    "project_to_los",
    "project_uncertainty_to_los",
    "read_epoch_displacements",
    "read_station_rate",
]
