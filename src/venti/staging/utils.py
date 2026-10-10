# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Utility functions for OPERA DISP-S1 processing.

This module provides common utility functions for date parsing, file parsing,
and database operations used across OPERA DISP-S1 processing scripts.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import requests  # type: ignore[import-untyped]
import rioxarray as rxr
from shapely.geometry import shape

logger = logging.getLogger(__name__)

# Date utilities
# Supported date formats for parsing
DATE_FORMATS = ["%Y-%m-%d", "%Y%m%d"]


def datetime_to_decimal_year(dt: datetime) -> float:
    """Convert datetime to decimal year with proper leap year handling.

    More accurate than 365.25-day approximation, especially for dates
    within a single year.

    Parameters
    ----------
    dt : datetime
        Calendar datetime to convert.

    Returns
    -------
    float
        Year expressed as decimal (e.g., 2014.5).

    Examples
    --------
    >>> datetime_to_decimal_year(datetime(2014, 1, 1))
    2014.0
    >>> datetime_to_decimal_year(datetime(2014, 7, 1))  # doctest: +SKIP
    2014.4945...

    """
    year = dt.year
    start_of_year = datetime(year, 1, 1)
    start_of_next_year = datetime(year + 1, 1, 1)

    year_elapsed = (dt - start_of_year).total_seconds()
    year_duration = (start_of_next_year - start_of_year).total_seconds()
    fraction = year_elapsed / year_duration

    return year + fraction


def parse_date(date_str: str | None) -> datetime | None:
    """Parse date string into datetime object.

    Supports ISO format (YYYY-MM-DD) and compact format (YYYYMMDD).

    Parameters
    ----------
    date_str : str or None
        Date string to parse, or None.

    Returns
    -------
    datetime or None
        Parsed datetime object, or None if input is None.

    Raises
    ------
    ValueError
        If date string format is invalid.

    Examples
    --------
    >>> parse_date("2024-01-15")
    datetime.datetime(2024, 1, 15, 0, 0)
    >>> parse_date("20240115")
    datetime.datetime(2024, 1, 15, 0, 0)
    >>> parse_date(None) is None
    True

    """
    if date_str is None:
        return None

    for fmt in DATE_FORMATS:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue

    msg = (
        f"Invalid date format: {date_str}. Expected formats: {', '.join(DATE_FORMATS)}"
    )
    raise ValueError(msg)


# Burst database utilities
# Default URL for OPERA consistent burst database
BURST_DB_URL = (
    "https://github.com/opera-adt/burst_db/releases/download/v0.13.0/"
    "opera-disp-s1-consistent-burst-ids-2025-09-16-2016-07-01_to_2024-12-31.json"
)


def load_burst_database(url: str = BURST_DB_URL) -> gpd.GeoDataFrame:
    """Load OPERA burst database JSON to GeoDataFrame.

    Parameters
    ----------
    url : str, optional
        URL to burst database JSON file.
        Default is OPERA consistent burst database.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with burst information including burst_id, geometry,
        and sensing times.

    Raises
    ------
    requests.HTTPError
        If database cannot be downloaded.

    Examples
    --------
    >>> gdf = load_burst_database()  # doctest: +SKIP
    >>> "burst_id" in gdf.columns  # doctest: +SKIP
    True

    """
    logger.info("Loading OPERA consistent burst database")

    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    logger.info(f"Database generation date: {data['metadata']['generation_time']}")

    features = []
    for burst_id, info in data["data"].items():
        features.append(
            {
                "frame_id": burst_id,
                "geometry": shape(info["geometry"]) if "geometry" in info else None,
                **{k: v for k, v in info.items() if k != "geometry"},
            }
        )

    return gpd.GeoDataFrame(features, crs="EPSG:4326")


# File parsing utilities
def combine_tropo_corrections(
    corrections_dir: Path,
    ref_time: datetime,
    sec_time: datetime,
) -> Path:
    """Combine per-epoch corrections into a single differential correction file.

    Subtracts the reference-date tropospheric delay from the secondary-date
    delay to form the net correction matching the sign convention of a DISP-S1
    displacement field (``secondary - reference``).

    Parameters
    ----------
    corrections_dir : Path
        Directory containing the per-epoch correction GeoTIFFs reprojected to
        the DISP UTM grid (e.g. ``tropo_corrections_32611/``).
    ref_time : datetime
        Reference epoch sensing time.
    sec_time : datetime
        Secondary epoch sensing time.

    Returns
    -------
    Path
        Path to the combined differential correction GeoTIFF, named
        ``tropo_correction_{ref_stem}_{sec_stem}_{epsg}.tif``.

    Raises
    ------
    FileNotFoundError
        If a correction file for either epoch cannot be found in
        `corrections_dir`.

    Examples
    --------
    >>> from datetime import datetime
    >>> combine_tropo_corrections(
    ...     Path("tropo/tropo_corrections_32611"),
    ...     datetime(2016, 7, 24, 1, 58, 9),
    ...     datetime(2016, 8, 5, 1, 58, 9),
    ... )  # doctest: +SKIP

    """
    ref_stem = ref_time.strftime("%Y%m%dT%H%M%S")
    sec_stem = sec_time.strftime("%Y%m%dT%H%M%S")

    ref_files = list(corrections_dir.glob(f"*{ref_stem}*.tif"))
    sec_files = list(corrections_dir.glob(f"*{sec_stem}*.tif"))

    if not ref_files:
        msg = (
            f"No correction file found for reference epoch {ref_stem} "
            f"in {corrections_dir}"
        )
        raise FileNotFoundError(msg)
    if not sec_files:
        msg = (
            f"No correction file found for secondary epoch {sec_stem} "
            f"in {corrections_dir}"
        )
        raise FileNotFoundError(msg)

    # EPSG suffix is the last underscore-delimited token in the directory name
    # e.g. "tropo_corrections_32611" -> "32611"
    epsg_suffix = corrections_dir.name.split("_")[-1]
    out_name = f"tropo_correction_{ref_stem}_{sec_stem}_{epsg_suffix}.tif"
    out_path = corrections_dir / out_name

    with (
        rxr.open_rasterio(ref_files[0]) as ref_da,
        rxr.open_rasterio(sec_files[0]) as sec_da,
    ):
        differential = sec_da - ref_da
        differential.rio.to_raster(out_path)

    logger.info(f"Differential tropo correction written to {out_path}")
    return out_path


def extract_sensing_times_from_file(disp_file: Path) -> list[datetime]:
    """Extract sensing times from DISP-S1 NetCDF filename.

    Parses DISP-S1 filename to extract reference and secondary dates.
    Expected format:
    OPERA_L3_DISP-S1_IW_F{frame}_VV_{ref_date}_{sec_date}_v{version}_{prod_date}.nc

    Parameters
    ----------
    disp_file : Path
        Path to DISP-S1 NetCDF file.

    Returns
    -------
    list[datetime]
        List of unique sensing times from reference and secondary dates.

    Raises
    ------
    FileNotFoundError
        If DISP file does not exist.
    ValueError
        If dates cannot be parsed from filename.

    Examples
    --------
    >>> from pathlib import Path
    >>> filename = Path("OPERA_L3_*.nc")
    >>> times = extract_sensing_times_from_file(filename)  # doctest: +SKIP
    >>> len(times)  # doctest: +SKIP
    2

    """
    if not disp_file.exists():
        msg = f"DISP file not found: {disp_file}"
        raise FileNotFoundError(msg)

    logger.info(f"Extracting sensing times from {disp_file.name}")

    filename = disp_file.name
    parts = filename.split("_")

    # Find parts matching datetime format (YYYYMMDDTHHMMSSZ)
    datetime_parts = [p for p in parts if len(p) == 16 and p.endswith("Z") and "T" in p]

    if len(datetime_parts) < 2:
        msg = (
            f"Cannot parse reference and secondary dates from filename: {filename}. "
            "Expected format: "
            "OPERA_L3_DISP-S1_IW_F{{frame}}_VV_{{ref_date}}_{{sec_date}}_"
            "v{{version}}_{{prod_date}}.nc"
        )
        raise ValueError(msg)

    ref_date_str = datetime_parts[0]
    sec_date_str = datetime_parts[1]

    try:
        ref_date = datetime.strptime(ref_date_str, "%Y%m%dT%H%M%SZ")
        sec_date = datetime.strptime(sec_date_str, "%Y%m%dT%H%M%SZ")
    except ValueError as e:
        msg = f"Failed to parse dates from filename {filename}: {e}"
        raise ValueError(msg) from e

    sensing_times = sorted({ref_date, sec_date})
    logger.info(
        f"Parsed sensing times: {ref_date.isoformat()} (ref), "
        f"{sec_date.isoformat()} (sec)"
    )

    return sensing_times


def extract_frame_id_from_filename(disp_file: Path) -> int:
    """Extract frame ID from DISP-S1 filename.

    Parses DISP-S1 filename to extract frame identifier.
    Expected format:
    OPERA_L3_DISP-S1_IW_F{frame}_VV_{ref_date}_{sec_date}_v{version}_{prod_date}.nc

    Parameters
    ----------
    disp_file : Path
        Path to DISP-S1 NetCDF file.

    Returns
    -------
    int
        Frame identifier.

    Raises
    ------
    ValueError
        If frame ID cannot be parsed from filename.

    Examples
    --------
    >>> from pathlib import Path
    >>> filename = Path("OPERA_L3_DISP-S1_*.nc")
    >>> extract_frame_id_from_filename(filename)  # doctest: +SKIP
    8887

    """
    filename = disp_file.name
    parts = filename.split("_")
    frame_parts = [p for p in parts if p.startswith("F") and p[1:].isdigit()]

    if not frame_parts:
        msg = (
            f"Cannot parse frame ID from filename: {filename}. "
            "Expected format: OPERA_L3_DISP-S1_IW_F{{frame}}_VV_..."
        )
        raise ValueError(msg)

    frame_str = frame_parts[0][1:]
    return int(frame_str)


def extract_sensing_times_from_database(
    frame_id: int,
    burst_db: gpd.GeoDataFrame,
    start: datetime | None = None,
    end: datetime | None = None,
) -> list[datetime]:
    """Extract sensing times for frame from burst database.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    burst_db : gpd.GeoDataFrame
        Burst database GeoDataFrame.
    start : datetime or None, optional
        Start date for filtering. Default is None (no start limit).
    end : datetime or None, optional
        End date for filtering. Default is None (no end limit).

    Returns
    -------
    list[datetime]
        List of unique sensing times for frame within date range.

    Raises
    ------
    ValueError
        If no bursts found for frame or no sensing times in date range.

    Examples
    --------
    >>> burst_db = load_burst_database()  # doctest: +SKIP
    >>> times = extract_sensing_times_from_database(8887, burst_db)  # doctest: +SKIP
    >>> len(times) > 0  # doctest: +SKIP
    True

    """
    logger.info(f"Extracting sensing times for frame {frame_id}")

    frame_bursts = burst_db[burst_db["frame_id"] == str(frame_id)]

    if frame_bursts.empty:
        msg = f"No bursts found for frame {frame_id} in database"
        raise ValueError(msg)

    sensing_times = []
    for _, burst in frame_bursts.iterrows():
        if "sensing_time_list" in burst:
            times = burst["sensing_time_list"]
            if times:
                sensing_times.extend(times)

    if not sensing_times:
        msg = f"No sensing times found for frame {frame_id} in database"
        raise ValueError(msg)

    # Convert to datetime
    unique_times = sorted({datetime.fromisoformat(t) for t in sensing_times})

    # Filter by date range
    if start:
        unique_times = [t for t in unique_times if t >= start]
    if end:
        unique_times = [t for t in unique_times if t <= end]

    if not unique_times:
        msg = (
            f"No sensing times found for frame {frame_id} "
            f"in date range [{start}, {end}]"
        )
        raise ValueError(msg)

    logger.info(f"Found {len(unique_times)} sensing times in date range")
    return unique_times
