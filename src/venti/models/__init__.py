# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
from .load_gia import (
    CARON_GIA,
    ICE6D_URL,
    clip_gia_df,
    download_ice6g_data,
    load_caron_model,
    load_ice6g_model,
    rasterize_gdf,
)

# Import ITRF data loading utilities
from .load_itrf import (
    ITRF14_DATA,
    ITRF20_DATA,
    convert_to_euler_poles,
    get_plate_data,
    json_to_dataframe,
    load_itrf_pmm,
)
from .plate_motion import (
    DEG_TO_RAD,
    # Constants
    EARTH_RADIUS_KM,
    MYR_TO_YEAR,
    RAD_TO_DEG,
    # Main calculation functions
    calculate_euler_pole,
    euler_pole_to_rotation_rate,
    # Coordinate system functions
    get_conversion_matrix,
    get_euler_pole_uncertainty,
    get_local_frame,
    # Prediction functions
    model_plate_velocities,
    model_velocities_from_euler_pole,
    # Conversion utilities
    rotation_rate_to_euler_pole,
)

# Package metadata
__all__ = [
    # Constants
    "CARON_GIA",
    "DEG_TO_RAD",
    "EARTH_RADIUS_KM",
    "ICE6D_URL",
    "ITRF14_DATA",
    "ITRF20_DATA",
    "MYR_TO_YEAR",
    "RAD_TO_DEG",
    # Euler pole calculation functions
    "calculate_euler_pole",
    # GIA clipping functions
    "clip_gia_df",
    # ITRF data loading functions
    "convert_to_euler_poles",
    # GIA data loading functions
    "download_ice6g_data",
    "euler_pole_to_rotation_rate",
    "get_conversion_matrix",
    "get_euler_pole_uncertainty",
    "get_local_frame",
    "get_plate_data",
    "json_to_dataframe",
    "load_caron_model",
    "load_ice6g_model",
    "load_itrf_pmm",
    # Modeling functions
    "model_plate_velocities",
    "model_velocities_from_euler_pole",
    "rasterize_gdf",
    "rotation_rate_to_euler_pole",
]
