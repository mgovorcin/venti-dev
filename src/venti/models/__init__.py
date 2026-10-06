# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
# Plate motion lives in geepers (single source, ADR-0010/ADR-0020); the
# names are re-exported here for the GIA/PMM context workflows.
from geepers.euler import (
    PLATE_CODES,
    EulerPole,
    load_plate_motion_model,
    plate_pole,
    plate_velocity_enu,
)

from .load_gia import (
    CARON_GIA,
    ICE6D_URL,
    clip_gia_df,
    download_ice6g_data,
    load_caron_model,
    load_ice6g_model,
    rasterize_gdf,
)

# Package metadata
__all__ = [
    # GIA
    "CARON_GIA",
    "ICE6D_URL",
    # plate motion (geepers)
    "PLATE_CODES",
    "EulerPole",
    "clip_gia_df",
    "download_ice6g_data",
    "load_caron_model",
    "load_ice6g_model",
    "load_plate_motion_model",
    "plate_pole",
    "plate_velocity_enu",
    "rasterize_gdf",
]
