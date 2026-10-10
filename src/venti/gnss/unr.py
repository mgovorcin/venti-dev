# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""UNR GNSS grid timeseries download and station velocity estimation.

Interfaces with the University of Nevada, Reno (UNR) geodesy database:
https://geodesy.unr.edu/grid_timeseries/
"""

from __future__ import annotations

import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests  # type: ignore[import-untyped]
from numpy.linalg import lstsq
from shapely.geometry import Point, box

logger = logging.getLogger(__name__)

# UNR timeseries endpoints keyed by reference frame.
#
# 'variable_data' holds per-epoch GNSS positions (time_variable_gridded); its
# sigma columns are position uncertainties (original GNSS uncertainty, plus
# gap-filling and spatial-interpolation-to-grid-point effects).
#
# 'constant_data' holds precomputed linear rates (time_contsant_gridded — sic,
# matches UNR's actual path spelling): east/north/up are exactly
# rate * (t - t0), and the sigma columns (constant across every row in a
# file) are the uncertainty of that rate, not of position. UNR only
# publishes this product for IGS20.
_UNR_URLS: dict[str, dict[str, str]] = {
    "IGS20": {
        "grid": (
            "https://geodesy.unr.edu/grid_timeseries/Version0.3/grid_latlon_lookup.txt"
        ),
        "variable_data": (
            "https://geodesy.unr.edu/grid_timeseries/Version0.3/time_variable_gridded/IGS20/"
        ),
        "constant_data": (
            "https://geodesy.unr.edu/grid_timeseries/Version0.3/time_contsant_gridded/IGS20/"
        ),
    },
    "IGS14": {
        "grid": "https://geodesy.unr.edu/grid_timeseries/grid_latlon_lookup.txt",
        "variable_data": (
            "https://geodesy.unr.edu/grid_timeseries/time_variable_gridded/IGS14/"
        ),
    },
}

# Column layout for UNR .tenv8 station files
_STATION_COLUMNS = ["year", "east", "north", "up", "sigma_e", "sigma_n", "sigma_u"]

# Epoch-displacement column names expected in GeoDataFrames returned by this module
_DISP_COLUMNS = ["deast", "dnorth", "dup", "dsigma_e", "dsigma_n", "dsigma_u"]


def download_grid_lookup(output_dir: Path, reference_frame: str = "IGS20") -> Path:
    """Download the UNR grid station lookup table.

    Parameters
    ----------
    output_dir : Path
        Directory to write the lookup file.
    reference_frame : str, optional
        GNSS reference frame, ``'IGS20'`` or ``'IGS14'``, by default ``'IGS20'``.

    Returns
    -------
    Path
        Path to the downloaded lookup file.

    Raises
    ------
    ValueError
        If `reference_frame` is not supported.
    RuntimeError
        If the download fails.

    Examples
    --------
    ::

        lookup = download_grid_lookup(Path("gnss"), reference_frame="IGS20")

    """
    if reference_frame not in _UNR_URLS:
        msg = (
            f"Unsupported reference frame '{reference_frame}'. Choose from"
            f" {list(_UNR_URLS)}"
        )
        raise ValueError(msg)

    url = _UNR_URLS[reference_frame]["grid"]
    output_path = output_dir / "grid_latlon_lookup.txt"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading UNR grid lookup from %s", url)
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    output_path.write_text(response.text, encoding="utf-8")

    logger.info("Saved grid lookup to %s", output_path)
    return output_path


def find_stations_in_bounds(
    grid_lookup_path: Path,
    bounds_snwe: tuple[float, float, float, float],
    utm_epsg: int,
    padding: float = 0.0,
) -> gpd.GeoDataFrame:
    """Filter the UNR grid station table to those within a bounding box.

    Parameters
    ----------
    grid_lookup_path : Path
        Path to the UNR ``grid_latlon_lookup.txt`` file.
    bounds_snwe : tuple of float
        Bounding box as ``(south, north, west, east)`` in the UTM CRS.
    utm_epsg : int
        EPSG code of the UTM projection matching `bounds_snwe`.
    padding : float, optional
        Extra padding in the same units as `bounds_snwe`, by default ``0.0``.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame of stations within the bounds, indexed by station ID,
        with columns ``lon``, ``lat``, and a ``geometry`` column in UTM.

    Examples
    --------
    ::

        stations = find_stations_in_bounds(
            lookup, bounds_snwe=(4.40e6, 4.64e6, 4.6e5, 7.5e5), utm_epsg=32618
        )

    """
    df = pd.read_csv(
        grid_lookup_path,
        sep=r"\s+",
        header=None,
        names=["id", "lon", "lat"],
    )
    gdf = gpd.GeoDataFrame(
        df,
        geometry=[
            Point(lon, lat) for lon, lat in zip(df["lon"], df["lat"], strict=False)
        ],
        crs="EPSG:4326",
    )
    gdf_utm = gdf.to_crs(epsg=utm_epsg)
    gdf_utm = gdf_utm.set_index("id")

    S, N, W, E = bounds_snwe
    bbox = box(W - padding, S - padding, E + padding, N + padding)
    return gdf_utm[gdf_utm.within(bbox)]


def download_station(
    station_id: int,
    output_dir: Path,
    reference_frame: str = "IGS20",
    grid_type: str = "constant",
) -> Path:
    """Download a single UNR GNSS grid station time series file.

    Parameters
    ----------
    station_id : int
        Numeric station ID (zero-padded to six digits in the filename).
    output_dir : Path
        Directory to write the station file.
    reference_frame : str, optional
        GNSS reference frame, by default ``'IGS20'``.
    grid_type : str, optional
        ``'constant'`` (precomputed linear rates, `time_contsant_gridded`) or
        ``'variable'`` (per-epoch positions, `time_variable_gridded`), by
        default ``'constant'``. UNR only publishes the constant-rate product
        for ``'IGS20'``.

    Returns
    -------
    Path
        Path to the downloaded station file,
        ``<id>_<reference_frame>_<grid_type>.tenv8``. Both products share one
        remote filename, so the grid type is kept in the local name to stop
        a cached file of one type being read as the other.

    Raises
    ------
    ValueError
        If `grid_type` is not available for `reference_frame`.
    RuntimeError
        If the download fails.

    Examples
    --------
    ::

        files = [
            download_station(sid, Path("gnss"), grid_type="constant")
            for sid in stations.index
        ]

    """
    filename = f"{station_id:06d}_{reference_frame}.tenv8"
    dest = output_dir / f"{station_id:06d}_{reference_frame}_{grid_type}.tenv8"

    if dest.exists():
        logger.debug("Station file already exists: %s", dest)
        return dest

    url_key = f"{grid_type}_data"
    if url_key not in _UNR_URLS[reference_frame]:
        msg = (
            f"grid_type='{grid_type}' is not available for reference_frame="
            f"'{reference_frame}'. UNR only publishes the constant-rate grid for IGS20."
        )
        raise ValueError(msg)

    url = _UNR_URLS[reference_frame][url_key] + filename
    logger.debug("Downloading station %s from %s", station_id, url)

    response = requests.get(url, timeout=60)
    if not response.ok:
        msg = f"Download failed for station {station_id} (HTTP {response.status_code})"
        raise RuntimeError(msg)

    dest.write_text(response.text, encoding="utf-8")
    return dest


def calculate_station_velocity(
    station_file: Path,
    start_year: float = 2014.0,
) -> tuple[float, float, float, float, float, float]:
    """Estimate east/north/up velocities via weighted linear regression.

    Fits a weighted least-squares model ``displacement = velocity * time + offset``
    to each component independently, using the inverse-variance weights.

    Parameters
    ----------
    station_file : Path
        UNR ``.tenv8`` station file (whitespace-delimited, columns ordered as
        ``year east north up sigma_e sigma_n sigma_u ...``).
    start_year : float, optional
        Exclude observations before this decimal year, by default ``2014.0``.

    Returns
    -------
    tuple of float
        ``(ve, vn, vu, sigma_ve, sigma_vn, sigma_vu)`` — velocities and their
        standard errors in the same units as the station file (mm/yr).

    Examples
    --------
    ::

        ve, vn, vu, se, sn, su = calculate_station_velocity(variable_file, 2016.0)

    """
    data = np.loadtxt(station_file)

    t = data[:, 0]
    east, north, up = data[:, 1], data[:, 2], data[:, 3]
    sigma_e, sigma_n, sigma_u = data[:, 4], data[:, 5], data[:, 6]

    mask = t >= start_year
    t, east, north, up = t[mask], east[mask], north[mask], up[mask]
    sigma_e, sigma_n, sigma_u = sigma_e[mask], sigma_n[mask], sigma_u[mask]

    # Design matrix [time, 1] for velocity + offset
    X = np.column_stack([t, np.ones_like(t)])

    def _wls(y: np.ndarray, sigma: np.ndarray) -> tuple[float, float]:
        """Weighted least-squares slope and its standard error."""
        w = 1.0 / sigma**2
        W = np.diag(w)
        XtW = X.T @ W
        beta, _, _, _ = lstsq(XtW @ X, XtW @ y, rcond=None)

        residuals = y - X @ beta
        dof = max(len(t) - 2, 1)
        rss = np.sum(w * residuals**2)
        cov = np.linalg.inv(XtW @ X)
        slope_std = np.sqrt(rss / dof * cov[0, 0])
        return float(beta[0]), float(slope_std)

    ve, ve_std = _wls(east, sigma_e)
    vn, vn_std = _wls(north, sigma_n)
    vu, vu_std = _wls(up, sigma_u)

    return ve, vn, vu, ve_std, vn_std, vu_std


def read_station_rate(
    station_file: Path,
) -> tuple[float, float, float, float, float, float]:
    """Read a precomputed east/north/up rate from a constant-grid station file.

    Constant-grid (`time_contsant_gridded`) files encode east/north/up as an
    exactly linear function of time, and the sigma columns as the
    uncertainty of that rate (constant across every row), rather than a
    per-epoch position uncertainty. The rate is therefore read directly
    rather than re-estimated: fitting residuals of an already-linear series
    would report near-zero uncertainty instead of UNR's propagated rate
    uncertainty.

    Parameters
    ----------
    station_file : Path
        UNR constant-grid ``.tenv8`` file (whitespace-delimited, columns
        ordered as ``year east north up sigma_e sigma_n sigma_u ...``).

    Returns
    -------
    tuple of float
        ``(ve, vn, vu, sigma_ve, sigma_vn, sigma_vu)`` — rates and their
        uncertainties in the same units as the station file (mm/yr).

    Examples
    --------
    ::

        ve, vn, vu, se, sn, su = read_station_rate(constant_file)  # mm/yr

    """
    data = np.loadtxt(station_file)

    t = data[:, 0]
    east, north, up = data[:, 1], data[:, 2], data[:, 3]
    sigma_ve, sigma_vn, sigma_vu = data[0, 4], data[0, 5], data[0, 6]

    ve = float(np.polyfit(t, east, 1)[0])
    vn = float(np.polyfit(t, north, 1)[0])
    vu = float(np.polyfit(t, up, 1)[0])

    return ve, vn, vu, float(sigma_ve), float(sigma_vn), float(sigma_vu)


def read_epoch_displacements(
    station_files: list[Path],
    ref_date: float,
    sec_date: float,
    station_gdf: gpd.GeoDataFrame,
    max_offset_days: float = 1.5,
) -> gpd.GeoDataFrame:
    """Read GNSS displacements for a specific epoch pair from station files.

    Uses the observation nearest to each date and returns ``sec - ref`` in
    east/north/up, the same sign as OPERA displacement and as the
    ``'constant'`` path (``rate * (sec_date - ref_date)``).

    Parameters
    ----------
    station_files : list of Path
        UNR station files to read.
    ref_date : float
        Reference epoch as decimal year.
    sec_date : float
        Secondary epoch as decimal year.
    station_gdf : gpd.GeoDataFrame
        GeoDataFrame of station metadata (geometry in UTM), indexed by station ID.
    max_offset_days : float, optional
        Skip a station (with a warning) if its nearest observation to either
        date is further away than this, by default 1.5 (UNR grids are daily).

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with columns ``deast``, ``dnorth``, ``dup``,
        ``dsigma_e``, ``dsigma_n``, ``dsigma_u``, and a ``geometry`` column.

    Raises
    ------
    ValueError
        If no station has observations near both dates.

    Examples
    --------
    ::

        disp_gdf = read_epoch_displacements(files, 2020.0, 2020.5, stations)
        gnss_los_mm = project_to_los(los_e, los_n, los_u, "product.nc", disp_gdf)

    """
    max_offset_years = max_offset_days / 365.25
    rows = []
    for path in station_files:
        station_id = int(path.name.split("_")[0])
        df = pd.read_csv(path, sep=r"\s+", names=_STATION_COLUMNS, usecols=range(7))

        ref_offset = (df["year"] - ref_date).abs()
        sec_offset = (df["year"] - sec_date).abs()
        if ref_offset.min() > max_offset_years or sec_offset.min() > max_offset_years:
            logger.warning(
                "Skipping GNSS station %d: no observation within %.1f days of "
                "%.4f and %.4f (record %.4f-%.4f)",
                station_id,
                max_offset_days,
                ref_date,
                sec_date,
                df["year"].iloc[0],
                df["year"].iloc[-1],
            )
            continue
        ref_row = df.iloc[ref_offset.argmin()]
        sec_row = df.iloc[sec_offset.argmin()]

        rows.append(
            {
                "id": station_id,
                "deast": sec_row["east"] - ref_row["east"],
                "dnorth": sec_row["north"] - ref_row["north"],
                "dup": sec_row["up"] - ref_row["up"],
                "dsigma_e": np.hypot(ref_row["sigma_e"], sec_row["sigma_e"]),
                "dsigma_n": np.hypot(ref_row["sigma_n"], sec_row["sigma_n"]),
                "dsigma_u": np.hypot(ref_row["sigma_u"], sec_row["sigma_u"]),
            }
        )

    if not rows:
        msg = f"No GNSS station covers both {ref_date:.4f} and {sec_date:.4f}"
        raise ValueError(msg)
    diff_df = pd.DataFrame(rows).set_index("id")
    return diff_df.join(station_gdf[["geometry"]], how="inner")
