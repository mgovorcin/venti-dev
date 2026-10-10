# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Data staging workflow for a single OPERA DISP-S1 frame.

Prepares data needed for a single DISP-S1 product:
  1. Download DISP-S1 product for the given frame and secondary date
  2. Generate GLO30 DEM for the frame extent
  3. Generate LOS ENU and incidence angle rasters
  4. Download and apply tropospheric corrections
  5. Download UNR GNSS station timeseries and compute velocities

Output directory structure::

    <output_dir>/
    ├── disp_s1/                    # Downloaded DISP-S1 NetCDF file
    ├── dem/                        # DEM in WGS84 and native UTM
    ├── los/                        # LOS ENU and incidence angle rasters
    ├── tropo/                      # Tropospheric correction products
    └── gnss/
        ├── grid_latlon_lookup.txt  # UNR station lookup table
        ├── stations/               # Per-station .tenv8 timeseries files
        └── velocities.parquet      # Station velocities (E/N/U, mm/yr)
"""

from __future__ import annotations

import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from opera_utils import get_frame_bbox

from venti.gnss.unr import (
    calculate_station_velocity,
    download_grid_lookup,
    download_station,
    find_stations_in_bounds,
    read_station_rate,
)

logger = logging.getLogger(__name__)

# Staging CLI scripts live outside the package; locate them relative to the repo root.
_STAGING_DIR = (
    Path(__file__).resolve().parent.parent.parent.parent / "scripts" / "staging"
)


def _find_disp_file(disp_dir: Path, frame_id: int) -> Path:
    """Locate the downloaded DISP-S1 file for a frame.

    Parameters
    ----------
    disp_dir : Path
        Directory containing downloaded DISP-S1 NetCDF files.
    frame_id : int
        OPERA frame identifier used to match filename.

    Returns
    -------
    Path
        Path to the matching DISP-S1 NetCDF file.

    Raises
    ------
    FileNotFoundError
        If no matching DISP-S1 file is found in ``disp_dir``.
    ValueError
        If multiple matching files are found and cannot be disambiguated.

    """
    pattern = f"OPERA_L3_DISP-S1*F{frame_id:05d}*.nc"
    matches = sorted(disp_dir.glob(pattern))

    if not matches:
        msg = (
            f"No DISP-S1 file found for frame {frame_id} in {disp_dir}. "
            f"Expected pattern: {pattern}"
        )
        raise FileNotFoundError(msg)

    if len(matches) > 1:
        msg = (
            f"Multiple DISP-S1 files found for frame {frame_id} in {disp_dir}: "
            f"{[m.name for m in matches]}. "
            "Narrow the date range to stage a single product."
        )
        raise ValueError(msg)

    return matches[0]


def download_gnss_data(
    frame_id: int,
    output_dir: Path,
    reference_frame: str = "IGS20",
    grid_type: str = "variable",
    padding: float = 0.0,
    num_workers: int = 4,
    start_year: float = 2014.0,
) -> Path:
    """Download UNR GNSS station timeseries and compute velocities for a frame.

    Finds all UNR grid stations within the frame's bounding box, downloads
    their ``.tenv8`` timeseries files in parallel, then estimates E/N/U
    velocities for each station and writes a parquet summary.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Root GNSS output directory. Lookup and velocity files are written here;
        station files go into ``output_dir/stations/``.
    reference_frame : str, optional
        GNSS reference frame (``'IGS20'`` or ``'IGS14'``). Default is ``'IGS20'``.
    grid_type : str, optional
        UNR grid product to download: ``'variable'`` (per-epoch positions,
        default — velocities are fit locally below via
        :func:`~venti.gnss.unr.calculate_station_velocity` over `start_year`
        onward) or ``'constant'`` (precomputed linear rates, IGS20 only —
        velocities are read directly via
        :func:`~venti.gnss.unr.read_station_rate`, and `start_year` does not
        apply).
    padding : float, optional
        Extra padding in meters beyond the frame extent when searching for
        stations. Default is ``0.0``.
    num_workers : int, optional
        Number of parallel download threads. Default is 4.
    start_year : float, optional
        Exclude observations before this decimal year when estimating
        velocities. Default is ``2014.0``. Only applies when
        ``grid_type='variable'``.

    Returns
    -------
    Path
        Path to the ``velocities.parquet`` file.

    Raises
    ------
    RuntimeError
        If velocity estimation fails for all stations.
    ValueError
        If no stations are found within the frame bounds, or if `grid_type`
        is not available for `reference_frame`.

    Examples
    --------
    Downloads from UNR, so not run as a doctest::

        vel_path = download_gnss_data(frame_id=8887, output_dir=Path("gnss"))
        vel_path = download_gnss_data(
            frame_id=8887, output_dir=Path("gnss"), grid_type="constant"
        )

    """
    output_dir.mkdir(parents=True, exist_ok=True)
    stations_dir = output_dir / "stations"
    stations_dir.mkdir(exist_ok=True)

    epsg, bbox = get_frame_bbox(frame_id)
    bounds_snwe = (bbox.bottom, bbox.top, bbox.left, bbox.right)

    lookup_path = download_grid_lookup(output_dir, reference_frame=reference_frame)

    station_gdf = find_stations_in_bounds(
        lookup_path, bounds_snwe, utm_epsg=epsg, padding=padding
    )
    if station_gdf.empty:
        msg = (
            f"No UNR stations found within frame {frame_id} bounds "
            f"(EPSG:{epsg}, padding={padding} m). "
            "Try increasing padding."
        )
        raise ValueError(msg)

    logger.info(
        "Found %d UNR stations within frame %d bounds", len(station_gdf), frame_id
    )

    def _download(sid: int) -> tuple[int, Path | None]:
        try:
            path = download_station(
                sid, stations_dir, reference_frame=reference_frame, grid_type=grid_type
            )
        except RuntimeError:
            logger.warning("Download failed for station %s — skipping", sid)
            return sid, None
        else:
            return sid, path

    station_files: dict[int, Path] = {}
    with ThreadPoolExecutor(max_workers=num_workers) as pool:
        futures = {pool.submit(_download, sid): sid for sid in station_gdf.index}
        for future in as_completed(futures):
            sid, path = future.result()
            if path is not None:
                station_files[sid] = path

    logger.info(
        "Downloaded %d / %d station files", len(station_files), len(station_gdf)
    )

    rows = []
    for sid, path in station_files.items():
        try:
            if grid_type == "constant":
                ve, vn, vu, ve_std, vn_std, vu_std = read_station_rate(path)
            else:
                ve, vn, vu, ve_std, vn_std, vu_std = calculate_station_velocity(
                    path, start_year=start_year
                )
            rows.append(
                {
                    "station_id": sid,
                    "lon": station_gdf.loc[sid, "lon"],
                    "lat": station_gdf.loc[sid, "lat"],
                    "geometry": station_gdf.loc[sid, "geometry"],
                    "ve_mmyr": ve,
                    "vn_mmyr": vn,
                    "vu_mmyr": vu,
                    "sigma_ve": ve_std,
                    "sigma_vn": vn_std,
                    "sigma_vu": vu_std,
                }
            )
        except Exception:
            logger.warning("Velocity estimation failed for station %s — skipping", sid)

    if not rows:
        msg = "Velocity estimation failed for all stations."
        raise RuntimeError(msg)

    import geopandas as gpd

    velocities_path = output_dir / "velocities.parquet"
    gpd.GeoDataFrame(rows, geometry="geometry").to_parquet(velocities_path, index=False)
    logger.info("Velocities for %d stations saved to %s", len(rows), velocities_path)

    return velocities_path


def stage_frame(
    frame_id: int,
    date: str,
    output_dir: Path,
    num_workers: int = 4,
    dem_buffer: float = 10_000.0,
    skip_tropo: bool = False,
    skip_gnss: bool = False,
    gnss_reference_frame: str = "IGS20",
    gnss_padding: float = 0.0,
    gnss_start_year: float = 2014.0,
) -> None:
    """Stage all ancillary data for a single OPERA DISP-S1 frame.

    Downloads the DISP-S1 product for the given secondary date and generates
    all ancillary data needed for downstream calibration: DEM, LOS geometry,
    tropospheric corrections, and UNR GNSS velocities.

    DEM and LOS generation are idempotent — existing outputs are reused if
    present.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    date : str
        Secondary date of the interferogram (YYYY-MM-DD or YYYYMMDD).
    output_dir : Path
        Root directory for all staged outputs.
    num_workers : int, optional
        Number of parallel workers for downloads and processing. Default is 4.
    dem_buffer : float, optional
        Buffer in meters around the frame extent for DEM generation.
        Default is 10,000 m (10 km).
    skip_tropo : bool, optional
        Skip tropospheric correction processing. Default is False.
    skip_gnss : bool, optional
        Skip UNR GNSS download and velocity estimation. Default is False.
    gnss_reference_frame : str, optional
        GNSS reference frame for UNR data (``'IGS20'`` or ``'IGS14'``).
        Default is ``'IGS20'``.
    gnss_padding : float, optional
        Extra padding in meters beyond frame bounds when searching for GNSS
        stations. Default is ``0.0``.
    gnss_start_year : float, optional
        Exclude GNSS observations before this decimal year when estimating
        velocities. Default is ``2014.0``.

    Raises
    ------
    FileNotFoundError
        If the downloaded DISP-S1 file cannot be located.
    ValueError
        If the frame ID is invalid or multiple products match the date.

    Examples
    --------
    Stage a single frame with all ancillary data::

        stage_frame(frame_id=8887, date="2016-06-15", output_dir=Path("./data"))

    Skip tropospheric corrections::

        stage_frame(
            frame_id=8887,
            date="2016-06-15",
            output_dir=Path("./data"),
            skip_tropo=True,
        )

    """
    sys.path.insert(0, str(_STAGING_DIR))
    from dem_cli import generate_frame_dem
    from disp_cli import download_frame_products
    from los_cli import generate_incidence_angle_raster, generate_los_enu_raster
    from tropo_cli import process_tropo_from_file
    from utils import parse_date

    sec_date = parse_date(date)
    assert sec_date is not None

    disp_dir = output_dir / "disp_s1"
    dem_dir = output_dir / "dem"
    los_dir = output_dir / "los"
    tropo_dir = output_dir / "tropo"
    gnss_dir = output_dir / "gnss"

    n_steps = 5 - int(skip_tropo) - int(skip_gnss)
    step = 0

    def _step(label: str) -> str:
        nonlocal step
        step += 1
        return f"[{step}/{n_steps}] {label}"

    logger.info(_step(f"Downloading DISP-S1 for frame {frame_id} on {date}"))
    download_frame_products(
        frame_id=frame_id,
        output_dir=disp_dir,
        start=sec_date,
        end=sec_date,
        num_workers=1,
    )
    disp_file = _find_disp_file(disp_dir, frame_id)
    logger.info("      Product: %s", disp_file.name)

    logger.info(_step(f"Generating DEM for frame {frame_id}"))
    dem_path = generate_frame_dem(
        frame_id=frame_id,
        buffer=dem_buffer,
        output_dir=dem_dir,
        use_disp_epsg=True,
    )
    logger.info("      DEM: %s", dem_path)

    logger.info(_step(f"Generating LOS geometry for frame {frame_id}"))
    los_path = generate_los_enu_raster(frame_id=frame_id, output_dir=los_dir)
    inc_path = generate_incidence_angle_raster(
        los_enu_path=los_path, output_dir=los_dir
    )
    logger.info("      LOS ENU: %s", los_path)
    logger.info("      Incidence angle: %s", inc_path)

    if not skip_tropo:
        logger.info(_step(f"Processing tropospheric corrections for {disp_file.name}"))
        process_tropo_from_file(
            disp_file=disp_file,
            dem_path=dem_path,
            incidence_angle_path=inc_path,
            output_dir=tropo_dir,
            num_workers=num_workers,
            to_disp_epsg=True,
        )
        logger.info("      Tropo corrections: %s", tropo_dir)

    if not skip_gnss:
        logger.info(_step(f"Downloading UNR GNSS data for frame {frame_id}"))
        vel_path = download_gnss_data(
            frame_id=frame_id,
            output_dir=gnss_dir,
            reference_frame=gnss_reference_frame,
            padding=gnss_padding,
            num_workers=num_workers,
            start_year=gnss_start_year,
        )
        logger.info("      Velocities: %s", vel_path)

    logger.info("Staging complete. All outputs in %s", output_dir.resolve())


def stage_window(
    frame_id: int,
    start: str,
    end: str,
    output_dir: Path,
    num_workers: int = 4,
    dem_buffer: float = 10_000.0,
    skip_tropo: bool = False,
    skip_gnss: bool = False,
    gnss_reference_frame: str = "IGS20",
    gnss_padding: float = 0.0,
    gnss_start_year: float = 2014.0,
) -> list[Path]:
    """Stage all ancillary data for multiple OPERA DISP-S1 products in a time window.

    Downloads every DISP-S1 product whose secondary date falls within
    ``[start, end]`` and generates all ancillary data required for downstream
    calibration.  Frame-level assets (DEM, LOS geometry, GNSS velocities) are
    produced once and shared across all products in the window.  Tropospheric
    corrections are downloaded in a single batched pass across all unique epoch
    sensing times and then combined into per-product differential correction
    files (secondary - reference).

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    start : str
        Start of the secondary-date window (YYYY-MM-DD or YYYYMMDD).
    end : str
        End of the secondary-date window (YYYY-MM-DD or YYYYMMDD).
    output_dir : Path
        Root directory for all staged outputs.  The directory layout mirrors
        that of :func:`stage_frame`:

        .. code-block:: text

            <output_dir>/
            ├── disp_s1/   # All downloaded DISP-S1 NetCDF files
            ├── dem/        # Shared GLO-30 DEM
            ├── los/        # Shared LOS ENU and incidence angle rasters
            ├── tropo/      # Per-product differential tropo corrections
            └── gnss/       # Shared UNR GNSS velocities

    num_workers : int, optional
        Number of parallel workers for downloads and processing. Default is 4.
    dem_buffer : float, optional
        Buffer in metres around the frame extent for DEM generation.
        Default is 10,000 m (10 km).
    skip_tropo : bool, optional
        Skip tropospheric correction processing. Default is False.
    skip_gnss : bool, optional
        Skip UNR GNSS download and velocity estimation. Default is False.
    gnss_reference_frame : str, optional
        GNSS reference frame for UNR data (``'IGS20'`` or ``'IGS14'``).
        Default is ``'IGS20'``.
    gnss_padding : float, optional
        Extra padding in metres beyond frame bounds when searching for GNSS
        stations. Default is ``0.0``.
    gnss_start_year : float, optional
        Exclude GNSS observations before this decimal year when estimating
        velocities. Default is ``2014.0``.

    Returns
    -------
    list[Path]
        Paths of the downloaded DISP-S1 NetCDF files, sorted by secondary date.

    Raises
    ------
    FileNotFoundError
        If no DISP-S1 files are found after download.
    ValueError
        If the frame ID is invalid or the date range yields no products.

    Examples
    --------
    Stage all products for frame 8887 between 2016-06-01 and 2016-12-31::

        stage_window(
            frame_id=8887,
            start="2016-06-01",
            end="2016-12-31",
            output_dir=Path("./data"),
        )

    Skip tropospheric corrections::

        stage_window(
            frame_id=8887,
            start="2016-06-01",
            end="2016-12-31",
            output_dir=Path("./data"),
            skip_tropo=True,
        )

    """
    sys.path.insert(0, str(_STAGING_DIR))
    from dem_cli import generate_frame_dem
    from disp_cli import download_frame_products
    from los_cli import generate_incidence_angle_raster, generate_los_enu_raster
    from tropo_cli import process_tropo_from_stack
    from utils import parse_date

    start_dt = parse_date(start)
    end_dt = parse_date(end)
    assert start_dt is not None
    assert end_dt is not None

    disp_dir = output_dir / "disp_s1"
    dem_dir = output_dir / "dem"
    los_dir = output_dir / "los"
    tropo_dir = output_dir / "tropo"
    gnss_dir = output_dir / "gnss"

    n_steps = 5 - int(skip_tropo) - int(skip_gnss)
    step = 0

    def _step(label: str) -> str:
        nonlocal step
        step += 1
        return f"[{step}/{n_steps}] {label}"

    logger.info(
        _step(
            f"Downloading DISP-S1 products for frame {frame_id} from {start} to {end}"
        )
    )
    download_frame_products(
        frame_id=frame_id,
        output_dir=disp_dir,
        start=start_dt,
        end=end_dt,
        num_workers=num_workers,
    )

    pattern = f"OPERA_L3_DISP-S1*F{frame_id:05d}*.nc"
    disp_files = sorted(disp_dir.glob(pattern))
    if not disp_files:
        msg = f"No DISP-S1 files found for frame {frame_id} in {disp_dir}"
        raise FileNotFoundError(msg)
    logger.info("      Found %d products", len(disp_files))

    logger.info(_step(f"Generating DEM for frame {frame_id}"))
    dem_path = generate_frame_dem(
        frame_id=frame_id,
        buffer=dem_buffer,
        output_dir=dem_dir,
        use_disp_epsg=True,
    )
    logger.info("      DEM: %s", dem_path)

    logger.info(_step(f"Generating LOS geometry for frame {frame_id}"))
    los_path = generate_los_enu_raster(frame_id=frame_id, output_dir=los_dir)
    inc_path = generate_incidence_angle_raster(
        los_enu_path=los_path, output_dir=los_dir
    )
    logger.info("      LOS ENU: %s", los_path)
    logger.info("      Incidence angle: %s", inc_path)

    if not skip_tropo:
        logger.info(
            _step(f"Processing tropospheric corrections for {len(disp_files)} products")
        )
        process_tropo_from_stack(
            disp_dir=disp_dir,
            dem_path=dem_path,
            incidence_angle_path=inc_path,
            output_dir=tropo_dir,
            num_workers=num_workers,
            frame_id=frame_id,
            to_disp_epsg=True,
        )
        logger.info("      Tropo corrections: %s", tropo_dir)

    if not skip_gnss:
        logger.info(_step(f"Downloading UNR GNSS data for frame {frame_id}"))
        vel_path = download_gnss_data(
            frame_id=frame_id,
            output_dir=gnss_dir,
            reference_frame=gnss_reference_frame,
            padding=gnss_padding,
            num_workers=num_workers,
            start_year=gnss_start_year,
        )
        logger.info("      Velocities: %s", vel_path)

    logger.info(
        "Staging complete. %d products staged in %s",
        len(disp_files),
        output_dir.resolve(),
    )
    return disp_files
