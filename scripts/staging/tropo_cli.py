# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tropospheric correction download and processing for DISP-S1 products.

This module provides tools to download, crop, and apply tropospheric corrections
from HRRR weather model data to DISP-S1 displacement products.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import asf_search as asf
import geopandas as gpd
import numpy as np
import requests  # type: ignore[import-untyped]
import rioxarray as rxr
import xarray as xr
from opera_utils import get_frame_bbox, get_frame_geojson
from opera_utils.tropo import apply_tropo, crop_tropo
from shapely.geometry import shape
from tqdm import tqdm
from utils import (
    extract_frame_id_from_filename,
    extract_sensing_times_from_database,
    extract_sensing_times_from_file,
    load_burst_database,
    parse_date,
)

# Constants
ASF_TROPO_COLLECTION = "C3717139408-ASF"
BURST_DB_URL = (
    "https://github.com/opera-adt/burst_db/releases/download/v0.13.0/"
    "opera-disp-s1-consistent-burst-ids-2025-09-16-2016-07-01_to_2024-12-31.json"
)
DATE_FORMATS = ("%Y-%m-%d", "%Y%m%d")
DEFAULT_OUTPUT_DIR = Path("./tropo")

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

# Suppress verbose logging from third-party libraries
logging.getLogger("asf_search").setLevel(logging.ERROR)
logging.getLogger("urllib3").setLevel(logging.WARNING)


def download_tropo_urls(output_dir: Path) -> Path:
    """Search and download list of tropospheric correction URLs.

    Parameters
    ----------
    output_dir : Path
        Directory to save the URL list.

    Returns
    -------
    Path
        Path to the saved tropo_urls.txt file.

    Examples
    --------
    >>> url_file = download_tropo_urls(Path("./tropo"))
    >>> print(url_file)
    ./tropo/tropo_urls.txt

    """
    logger.info("Searching ASF for tropospheric corrections")

    results = asf.search(collections=[ASF_TROPO_COLLECTION])

    if not results:
        msg = (
            "No tropospheric correction products found in collection"
            f" {ASF_TROPO_COLLECTION}"
        )
        raise ValueError(msg)

    logger.info(f"Found {len(results)} tropospheric correction products")

    tropo_gdf = gpd.GeoDataFrame.from_features(results.geojson())

    if "url" not in tropo_gdf.columns:
        msg = "URL column not found in search results"
        raise ValueError(msg)

    output_dir.mkdir(parents=True, exist_ok=True)
    url_file = output_dir / "tropo_urls.txt"

    tropo_gdf["url"].to_csv(url_file, header=False, index=False)
    logger.info(f"Saved {len(tropo_gdf)} URLs to {url_file}")

    return url_file


def download_prep_tropo(
    sensing_times: list[datetime],
    frame_bounds: tuple[float, float, float, float],
    dem_path: Path,
    incidence_angle_path: Path,
    output_dir: Path,
    num_workers: int = 2,
    dem_hgt_buffer: float = 3000.0,
) -> None:
    """Download and prepare tropospheric correction data.

    Parameters
    ----------
    sensing_times : list[datetime]
        List of sensing times to process.
    frame_bounds : tuple[float, float, float, float]
        Geographic bounds (minx, miny, maxx, maxy).
    dem_path : Path
        Path to DEM file.
    incidence_angle_path : Path
        Path to incidence angle file.
    output_dir : Path
        Output directory for all products.
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.

    Raises
    ------
    FileNotFoundError
        If DEM or incidence angle file does not exist.

    """
    if not dem_path.exists():
        msg = f"DEM file not found: {dem_path}"
        raise FileNotFoundError(msg)

    if not incidence_angle_path.exists():
        msg = f"Incidence angle file not found: {incidence_angle_path}"
        raise FileNotFoundError(msg)

    # Download URL list
    url_file = download_tropo_urls(output_dir)

    # Crop to frame and times
    logger.info(f"Cropping tropospheric data for {len(sensing_times)} sensing times")
    cropped_dir = output_dir / "cropped_tropo"
    cropped_dir.mkdir(parents=True, exist_ok=True)

    # Get max height from DEM with buffer for tropo download/crop
    with rxr.open_rasterio(dem_path) as hgt:
        max_dem_hgt = float(hgt.squeeze().max().values)
        hgt_max = max_dem_hgt + dem_hgt_buffer

    crop_tropo(
        url_file,
        sensing_times,
        frame_bounds,
        height_max=hgt_max,
        output_dir=cropped_dir,
        num_workers=num_workers,
    )
    cropped_files = sorted(cropped_dir.glob("tropo*.nc"))
    logger.info(f"Generated {len(cropped_files)} cropped files")

    # Apply corrections
    corrections_dir = output_dir / "tropo_corrections"
    corrections_dir.mkdir(parents=True, exist_ok=True)

    apply_tropo(
        cropped_files,
        dem_path=dem_path,
        incidence_angle_path=incidence_angle_path,
        output_dir=corrections_dir,
        subtract_first_date=False,  # Skip, do it later
        num_workers=num_workers,
    )

    logger.info(f"Corrections saved to {corrections_dir}")


def _reproject_to_disp_utm(input_file: Path, disp_file: Path) -> Path:
    """Reproject tropospheric correction to match DISP-S1 UTM grid.

    Parameters
    ----------
    input_file : Path
        Path to input tropospheric correction file.
    disp_file : Path
        Path to DISP-S1 file to match projection.

    Returns
    -------
    Path
        Path to reprojected output file.

    """
    with rxr.open_rasterio(input_file) as src, rxr.open_rasterio(disp_file) as disp:
        out = src.rio.reproject_match(disp)
        epsg = disp.rio.crs.to_epsg()
    out = out.squeeze(drop=True)
    output_dir = input_file.parents[1] / f"tropo_corrections_{epsg}"
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{input_file.stem}_{epsg}.tif"
    output_path = output_dir / filename

    out.rio.to_raster(output_path)

    return output_path


def _reproject_to_epsg(input_file: Path, epsg: int) -> Path:
    """Reproject tropospheric correction to specified EPSG code.

    Parameters
    ----------
    input_file : Path
        Path to input tropospheric correction file.
    epsg : int
        Target EPSG code.

    Returns
    -------
    Path
        Path to reprojected output file.

    """
    with rxr.open_rasterio(input_file) as src:
        out = src.rio.reproject(f"EPSG:{epsg}")

    output_dir = input_file.parents[1] / f"tropo_corrections_{epsg}"
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{input_file.stem}_{epsg}.tif"
    output_path = output_dir / filename

    out.rio.to_raster(output_path)

    return output_path


def _reproject_files_parallel(
    tropo_files: list[Path],
    reproject_fn: Callable,
    reproject_args: tuple,
    num_workers: int = 2,
) -> list[Path]:
    """Reproject multiple files in parallel.

    Parameters
    ----------
    tropo_files : list[Path]
        List of tropospheric correction files to reproject.
    reproject_fn : Callable
        Reprojection function to apply.
    reproject_args : tuple
        Additional arguments to pass to reprojection function.
    num_workers : int, optional
        Number of parallel workers. Default is 2.

    Returns
    -------
    list[Path]
        List of reprojected output file paths.

    """
    output_paths = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(reproject_fn, file, *reproject_args): file
            for file in tropo_files
        }

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Reprojecting files",
            unit="file",
        ):
            input_file = futures[future]
            try:
                output_path = future.result()
                output_paths.append(output_path)
            except Exception:
                logger.exception(f"Failed to reproject {input_file}")
                raise
    logger.info(f"Corrections saved to {Path(output_paths[0]).parent}")
    return sorted(output_paths)


def process_tropo_from_file(
    disp_file: Path,
    dem_path: Path,
    incidence_angle_path: Path,
    output_dir: Path,
    num_workers: int = 2,
    dem_hgt_buffer: float = 3000.0,
    to_disp_epsg: bool = True,
) -> None:
    """Process tropospheric corrections using a DISP-S1 file.

    Extracts sensing times from the DISP file, downloads and crops
    tropospheric data, then applies corrections.

    Parameters
    ----------
    disp_file : Path
        Path to DISP-S1 NetCDF file.
    dem_path : Path
        Path to DEM file.
    incidence_angle_path : Path
        Path to incidence angle file.
    output_dir : Path
        Output directory for all products.
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.
    to_disp_epsg : bool, optional
        Reproject to DISP file's UTM projection. Default is True.

    Examples
    --------
    >>> process_tropo_from_file(
    ...     disp_file=Path("DISP_frame_8887.nc"),
    ...     dem_path=Path("dem.tif"),
    ...     incidence_angle_path=Path("incidence.tif"),
    ...     output_dir=Path("./tropo")
    ... )

    """
    logger.info(f"Processing tropospheric corrections from file: {disp_file}")

    sensing_times = extract_sensing_times_from_file(disp_file)
    frame_id = extract_frame_id_from_filename(disp_file)

    selected_frame = get_frame_geojson([frame_id], as_geodataframe=True)
    frame_bounds = tuple(selected_frame.bounds.values[0])

    download_prep_tropo(
        sensing_times,
        frame_bounds,
        dem_path,
        incidence_angle_path,
        output_dir,
        num_workers,
        dem_hgt_buffer,
    )

    if to_disp_epsg:
        tropo_files = list((output_dir / "tropo_corrections").glob("*.tif"))
        logger.info(f"Reprojecting {len(tropo_files)} files to DISP grid in parallel")
        _reproject_files_parallel(
            tropo_files, _reproject_to_disp_utm, (disp_file,), num_workers
        )


def process_tropo_from_stack(
    disp_dir: Path,
    dem_path: Path,
    incidence_angle_path: Path,
    output_dir: Path,
    num_workers: int = 2,
    dem_hgt_buffer: float = 3000.0,
    frame_id: int | None = None,
    to_disp_epsg: bool = True,
) -> None:
    """Process tropospheric corrections from a directory of DISP-S1 files.

    Extracts sensing times from all DISP files in the directory,
    downloads and crops tropospheric data, then applies corrections.

    Parameters
    ----------
    disp_dir : Path
        Path to directory containing DISP-S1 NetCDF files.
    dem_path : Path
        Path to DEM file.
    incidence_angle_path : Path
        Path to incidence angle file.
    output_dir : Path
        Output directory for all products.
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.
    frame_id : int or None, optional
        Specific frame ID to process. If None, processes all frames.
        Default is None.
    to_disp_epsg : bool, optional
        Reproject to DISP file's UTM projection. Default is True.

    Raises
    ------
    ValueError
        If multiple frame IDs are detected without specifying which to process.

    Examples
    --------
    >>> process_tropo_from_stack(
    ...     disp_dir=Path("./disp_files"),
    ...     dem_path=Path("dem.tif"),
    ...     incidence_angle_path=Path("incidence.tif"),
    ...     output_dir=Path("./tropo"),
    ...     frame_id=8887
    ... )

    """
    logger.info(f"Processing tropospheric corrections from dir: {disp_dir}")
    disp_list = list(Path(disp_dir).glob("OPERA_L3_DISP-S1*.nc"))
    logger.info(f"Found {len(disp_list)} DISP-S1 files in {disp_dir}")

    frame_id_list = [extract_frame_id_from_filename(f) for f in disp_list]
    unique_frame_ids = np.unique(frame_id_list)

    # Filter by frame_id if specified
    if frame_id is not None:
        disp_list = [
            f
            for f, fid in zip(disp_list, frame_id_list, strict=True)
            if fid == frame_id
        ]
        target_frame_id = frame_id
    else:
        if len(unique_frame_ids) > 1:
            msg = (
                f"Multiple frame IDs detected in {disp_dir}: {unique_frame_ids}. "
                "Rerun with specified frame_id parameter."
            )
            raise ValueError(msg)
        target_frame_id = int(unique_frame_ids[0])

    # Extract and flatten sensing times
    all_sensing_times = []
    for file in disp_list:
        times = extract_sensing_times_from_file(file)
        all_sensing_times.append(times)
    sensing_times = np.unique(all_sensing_times).tolist()

    selected_frame = get_frame_geojson([target_frame_id], as_geodataframe=True)
    frame_bounds = tuple(selected_frame.bounds.values[0])

    download_prep_tropo(
        sensing_times,
        frame_bounds,
        dem_path,
        incidence_angle_path,
        output_dir,
        num_workers,
        dem_hgt_buffer,
    )

    if to_disp_epsg:
        tropo_files = list((output_dir / "tropo_corrections").glob("*.tif"))
        logger.info(f"Reprojecting {len(tropo_files)} files to DISP grid in parallel")
        _reproject_files_parallel(
            tropo_files, _reproject_to_disp_utm, (disp_list[0],), num_workers
        )


def process_tropo_from_db(
    frame_id: int,
    dem_path: Path,
    incidence_angle_path: Path,
    output_dir: Path,
    num_workers: int = 2,
    dem_hgt_buffer: float = 3000.0,
    start: datetime | None = None,
    end: datetime | None = None,
    to_disp_epsg: bool = True,
) -> None:
    """Process tropospheric corrections using frame ID and date range.

    Queries the OPERA consistent burst database for sensing times,
    downloads and crops tropospheric data, then applies corrections.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    dem_path : Path
        Path to DEM file.
    incidence_angle_path : Path
        Path to incidence angle file.
    output_dir : Path
        Output directory for all products.
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.
    start : datetime or None, optional
        Start date for filtering. Default is None (no start limit).
    end : datetime or None, optional
        End date for filtering. Default is None (no end limit).
    to_disp_epsg : bool, optional
        Reproject to frame's native EPSG. Default is True.

    Examples
    --------
    >>> process_tropo_from_db(
    ...     frame_id=8887,
    ...     dem_path=Path("dem.tif"),
    ...     incidence_angle_path=Path("incidence.tif"),
    ...     output_dir=Path("./tropo"),
    ...     start=datetime(2024, 1, 1),
    ...     end=datetime(2024, 12, 31)
    ... )

    """
    logger.info(
        f"Processing tropospheric corrections for frame {frame_id} "
        f"from {start or 'beginning'} to {end or 'present'}"
    )

    burst_db = load_burst_database()
    sensing_times = extract_sensing_times_from_database(frame_id, burst_db, start, end)

    selected_frame = get_frame_geojson([frame_id], as_geodataframe=True)
    frame_bounds = tuple(selected_frame.bounds.values[0])
    epsg = get_frame_bbox(frame_id)[0]

    download_prep_tropo(
        sensing_times,
        frame_bounds,
        dem_path,
        incidence_angle_path,
        output_dir,
        num_workers,
        dem_hgt_buffer,
    )

    if to_disp_epsg:
        tropo_files = list((output_dir / "tropo_corrections").glob("*.tif"))
        logger.info(f"Reprojecting {len(tropo_files)} files to EPSG:{epsg} in parallel")
        _reproject_files_parallel(tropo_files, _reproject_to_epsg, (epsg,), num_workers)


@dataclass
class FromFile:
    """Process tropospheric corrections from a DISP-S1 file.

    Attributes
    ----------
    disp_file : Path
        Path to DISP-S1 NetCDF file containing reference and secondary dates.
    dem_path : Path
        Path to DEM file for the frame.
    incidence_angle_path : Path
        Path to incidence angle file for the frame.
    output_dir : Path
        Output directory for tropospheric products.
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.
    to_disp_epsg : bool, optional
        Reproject to DISP file's UTM projection. Default is True.

    """

    disp_file: Path
    dem_path: Path
    incidence_angle_path: Path
    output_dir: Path = DEFAULT_OUTPUT_DIR
    num_workers: int = 2
    dem_hgt_buffer: float = 3000.0
    to_disp_epsg: bool = True

    def __call__(self) -> None:
        """Execute file-based tropospheric correction processing."""
        process_tropo_from_file(
            self.disp_file,
            self.dem_path,
            self.incidence_angle_path,
            self.output_dir,
            self.num_workers,
            self.dem_hgt_buffer,
            self.to_disp_epsg,
        )


@dataclass
class FromStack:
    """Process tropospheric corrections from directory with DISP-S1 files.

    Attributes
    ----------
    disp_dir : Path
        Path to directory containing DISP-S1 NetCDF files.
    dem_path : Path
        Path to DEM file for the frame.
    incidence_angle_path : Path
        Path to incidence angle file for the frame.
    output_dir : Path
        Output directory for tropospheric products.
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.
    frame_id : int or None, optional
        Specific frame ID to process. If None, processes all frames.
        Default is None.
    to_disp_epsg : bool, optional
        Reproject to DISP file's UTM projection. Default is True.

    """

    disp_dir: Path
    dem_path: Path
    incidence_angle_path: Path
    output_dir: Path = DEFAULT_OUTPUT_DIR
    num_workers: int = 2
    dem_hgt_buffer: float = 3000.0
    frame_id: int | None = None
    to_disp_epsg: bool = True

    def __call__(self) -> None:
        """Execute directory-based tropospheric correction processing."""
        process_tropo_from_stack(
            self.disp_dir,
            self.dem_path,
            self.incidence_angle_path,
            self.output_dir,
            self.num_workers,
            self.dem_hgt_buffer,
            self.frame_id,
            self.to_disp_epsg,
        )


@dataclass
class FromDatabase:
    """Process tropospheric corrections from frame ID and date range.

    Attributes
    ----------
    frame_id : int
        OPERA frame identifier.
    dem_path : Path
        Path to DEM file for the frame.
    incidence_angle_path : Path
        Path to incidence angle file for the frame.
    output_dir : Path
        Output directory for tropospheric products.
    start : str or None, optional
        Start date (YYYY-MM-DD or YYYYMMDD format).
    end : str or None, optional
        End date (YYYY-MM-DD or YYYYMMDD format).
    num_workers : int, optional
        Number of parallel workers. Default is 2.
    dem_hgt_buffer : float, optional
        Height buffer above DEM max (meters). Default is 3000.0.
    to_disp_epsg : bool, optional
        Reproject to frame's native EPSG. Default is True.

    """

    frame_id: int
    dem_path: Path
    incidence_angle_path: Path
    output_dir: Path = DEFAULT_OUTPUT_DIR
    start: str | None = None
    end: str | None = None
    num_workers: int = 2
    dem_hgt_buffer: float = 3000.0
    to_disp_epsg: bool = True

    def __call__(self) -> None:
        """Execute database-based tropospheric correction processing."""
        process_tropo_from_db(
            self.frame_id,
            self.dem_path,
            self.incidence_angle_path,
            self.output_dir,
            num_workers=self.num_workers,
            dem_hgt_buffer=self.dem_hgt_buffer,
            start=parse_date(self.start),
            end=parse_date(self.end),
            to_disp_epsg=self.to_disp_epsg,
        )


def main() -> None:
    """Run CLI with file and stack subcommands."""
    import tyro

    tyro.extras.subcommand_cli_from_dict(
        {"db": FromDatabase, "file": FromFile, "stack": FromStack}
    )()


if __name__ == "__main__":
    main()
