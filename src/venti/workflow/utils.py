# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Utility functions for timeseries calibration.

This module provides helper functions for date handling, downsampling,
and other common operations.
"""

from __future__ import annotations

import functools
import logging
import os
import re
import shutil
import tempfile
from collections.abc import Callable, Iterator
from contextlib import contextmanager, suppress
from datetime import date, datetime
from pathlib import Path
from typing import Any, TypeVar

import numpy as np

logger = logging.getLogger(__name__)


def extract_dates_from_filename(
    filename: str | Path,
) -> tuple[date | None, date | None]:
    """Extract reference and secondary dates from OPERA filename.

    Parameters
    ----------
    filename : str or Path
        Filename containing dates in format YYYYMMDDTHHMMSS

    Returns
    -------
    tuple
        (reference_date, secondary_date) as datetime.date objects
        Returns (None, date) if only one date found

    Examples
    --------
    ::

        ref, sec = extract_dates_from_filename(
            'OPERA_L3_DISP-S1_20200101T000000_20200115T000000.nc'
        )

    """
    filename = str(filename)
    pattern = r"\d{8}T\d{6}"
    matches = re.findall(pattern, filename)

    start_date, end_date = None, None

    if len(matches) >= 2:
        start_date = datetime.strptime(matches[0], "%Y%m%dT%H%M%S").date()
        end_date = datetime.strptime(matches[1], "%Y%m%dT%H%M%S").date()
    elif len(matches) == 1:
        end_date = datetime.strptime(matches[0], "%Y%m%dT%H%M%S").date()

    return start_date, end_date


def datetime_to_decimal_year(dt: datetime) -> float:
    """Convert datetime to decimal year.

    Parameters
    ----------
    dt : datetime
        DateTime object to convert

    Returns
    -------
    float
        Decimal year (e.g., 2020.5 for middle of 2020)

    Examples
    --------
    ::

        dec_year = datetime_to_decimal_year(datetime(2020, 7, 1))
        # Returns approximately 2020.5

    """
    year_start = datetime(dt.year, 1, 1)
    next_year_start = datetime(dt.year + 1, 1, 1)
    year_length = (next_year_start - year_start).total_seconds()
    seconds_into_year = (dt - year_start).total_seconds()
    decimal_year = dt.year + seconds_into_year / year_length
    return round(decimal_year, 4)


def match_correction_to_displacement(
    correction_files: list[Path] | None, displacement_files: list[Path]
) -> list[tuple[Path | None, Path | None, Path]]:
    """Match per-epoch correction files to displacement files by date.

    Each displacement file covers a reference-secondary epoch pair.  The
    caller is expected to supply one correction file per sensing epoch so
    that the calibration workflow can form the differential correction
    (secondary minus reference) itself.

    Parameters
    ----------
    correction_files : list of Path or None
        Per-epoch correction files (one file per sensing time, containing a
        single date in the filename).  If ``None``, returns triples with
        ``None`` for both correction slots.
    displacement_files : list of Path
        Displacement files, each containing both a reference and a secondary
        date in the filename.

    Returns
    -------
    list of tuple
        ``(ref_correction_file, sec_correction_file, displacement_file)``
        triples.  Either correction slot is ``None`` when the corresponding
        epoch file cannot be found.

    Examples
    --------
    ::

        tropo_files = sorted(Path('tropo/tropo_corrections_32611/').glob('*.tif'))
        disp_files = sorted(Path('disp/').glob('*.nc'))
        matches = match_correction_to_displacement(tropo_files, disp_files)

    """
    disp_dict = {extract_dates_from_filename(f): f for f in displacement_files}

    if correction_files is None:
        return [(None, None, disp_file) for disp_file in disp_dict.values()]

    # Build a lookup keyed by the single epoch date found in each correction file
    corr_dict: dict[date, Path] = {}
    for f in correction_files:
        ref_date, sec_date = extract_dates_from_filename(f)
        epoch_date = sec_date if sec_date is not None else ref_date
        if epoch_date is not None:
            corr_dict[epoch_date] = f

    matches: list[tuple[Path | None, Path | None, Path]] = []
    for (ref_date, sec_date), disp_file in disp_dict.items():
        ref_tropo = corr_dict.get(ref_date) if ref_date is not None else None
        sec_tropo = corr_dict.get(sec_date) if sec_date is not None else None
        if ref_tropo is None or sec_tropo is None:
            logger.warning(
                f"Missing tropo file(s) for {disp_file.name}: "
                f"ref={'found' if ref_tropo else 'missing'}, "
                f"sec={'found' if sec_tropo else 'missing'}"
            )
        matches.append((ref_tropo, sec_tropo, disp_file))

    return matches


def compute_average_temporal_coherence(
    netcdf_files: list[Path], output_dir: Path, variable: str = "temporal_coherence"
) -> Path:
    """Compute average temporal coherence from multiple NetCDF files.

    Parameters
    ----------
    netcdf_files : list of Path
        List of NetCDF files
    output_dir : Path
        Output directory
    variable : str, optional
        Variable name, by default 'temporal_coherence'

    Returns
    -------
    Path
        Path to output GeoTIFF with average temporal coherence

    Examples
    --------
    ::

        nc_files = sorted(Path('data/').glob('*.nc'))
        avg_coh = compute_average_temporal_coherence(nc_files, Path('output/'))

    """
    import rioxarray  # noqa: F401  (registers the .rio accessor used below)
    import xarray as xr

    output_file = output_dir / f"average_{variable}.tif"

    if output_file.exists():
        logger.info(f"Using existing average {variable}: {output_file}")
        return output_file

    logger.info(f"Computing average {variable} from {len(netcdf_files)} files")

    # Grid coordinates and CRS come from the first product; without them the
    # GeoTIFF would be written with an identity (pixel-index) transform.
    crs = None
    with xr.open_dataset(netcdf_files[0]) as ds:
        coords = {"y": ds["y"].values, "x": ds["x"].values}
        if "spatial_ref" in ds:
            sr = ds["spatial_ref"]
            for attr in ("crs_wkt", "spatial_ref", "wkt"):
                if attr in sr.attrs:
                    crs = sr.attrs[attr]
                    break
        elif "crs" in ds.attrs:
            crs = ds.attrs["crs"]

    if crs is None:
        logger.warning("Could not extract CRS from NetCDF")

    # Compute incremental mean to avoid loading all files into memory at once
    acc = None
    count = 0
    for f in netcdf_files:
        with xr.open_dataset(f) as ds:
            if variable not in ds:
                continue
            data = ds[variable].values
        if acc is None:
            acc = data.astype("float64")
        else:
            acc += data
        count += 1

    if acc is None or count == 0:
        msg = f"Variable '{variable}' not found in any files"
        raise ValueError(msg)

    average_da = xr.DataArray(
        (acc / count).astype("float32"),
        dims=["y", "x"],
        coords=coords,
    )

    if crs is not None:
        average_da = average_da.rio.write_crs(crs)
        average_da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    average_da.rio.to_raster(output_file)

    logger.info(f"Saved average {variable} to: {output_file}")
    return output_file


def parse_window_size_meters(window_meters: float, posting_meters: float = 30.0) -> int:
    """Convert window size from meters to pixels.

    Parameters
    ----------
    window_meters : float
        Window size in meters
    posting_meters : float, optional
        Pixel posting in meters, by default 30.0 (for OPERA)

    Returns
    -------
    int
        Window size in pixels

    Examples
    --------
    ::

        win_size = parse_window_size_meters(30000, posting_meters=30)
        # Returns 1000 pixels

    """
    window_pixels = int(np.rint(window_meters / posting_meters))
    return window_pixels


def get_file_dates(file_path: str | Path) -> tuple[float, float]:
    """Extract dates from file and convert to decimal years.

    Parameters
    ----------
    file_path : str or Path
        Path to file with dates in filename

    Returns
    -------
    tuple
        (ref_decimal_year, sec_decimal_year)

    Examples
    --------
    ::

        ref_year, sec_year = get_file_dates('OPERA_20200101T000000_20200115T000000.nc')

    """
    from opera_utils import get_dates

    dates = get_dates(str(file_path))
    ref_decimal = datetime_to_decimal_year(dates[0])
    sec_decimal = datetime_to_decimal_year(dates[1])

    return ref_decimal, sec_decimal


@contextmanager
def scratch_temp_dir(scratch_path: str | Path) -> Iterator[Path]:
    """Send all temporary files to a private folder under ``<scratch_path>/tmp``.

    Covers Python's `tempfile` (Venti's own temporary files and any library
    using it) and GDAL work files (``CPL_TMPDIR``), so a run never depends
    on the working directory or ``/tmp`` being writable. Each call gets its
    own ``<scratch_path>/tmp/run_*`` folder, so concurrent runs sharing a
    scratch path never touch each other's files; the folder is deleted and
    previous settings are restored on exit.

    Parameters
    ----------
    scratch_path : str or Path
        Scratch directory from ``product_path_group.scratch_path``.

    Yields
    ------
    Path
        This run's temporary directory.

    Raises
    ------
    OSError
        If ``<scratch_path>/tmp`` cannot be created.

    Examples
    --------
    >>> import tempfile
    >>> from venti.workflow.utils import scratch_temp_dir
    >>> with scratch_temp_dir(tempfile.mkdtemp()) as tmp:
    ...     tempfile.gettempdir() == str(tmp) and tmp.parent.name == "tmp"
    True

    """
    tmp_root = Path(scratch_path).resolve() / "tmp"
    try:
        tmp_root.mkdir(parents=True, exist_ok=True)
        tmp_dir = Path(tempfile.mkdtemp(prefix="run_", dir=tmp_root))
    except OSError as exc:
        msg = (
            f"Cannot create a temporary directory in {tmp_root}; set "
            "product_path_group.scratch_path in the runconfig to a writable directory"
        )
        raise OSError(msg) from exc

    previous_tempdir = tempfile.tempdir
    previous_cpl = os.environ.get("CPL_TMPDIR")
    tempfile.tempdir = str(tmp_dir)
    os.environ["CPL_TMPDIR"] = str(tmp_dir)
    logger.info("Temporary files go to %s", tmp_dir)
    try:
        yield tmp_dir
    finally:
        tempfile.tempdir = previous_tempdir
        if previous_cpl is None:
            os.environ.pop("CPL_TMPDIR", None)
        else:
            os.environ["CPL_TMPDIR"] = previous_cpl
        shutil.rmtree(tmp_dir, ignore_errors=True)
        with suppress(OSError):
            tmp_root.rmdir()  # only succeeds once no other run is using it


_Method = TypeVar("_Method", bound=Callable[..., Any])


def with_scratch_temp_dir(method: _Method) -> _Method:
    """Run a workflow method inside `scratch_temp_dir` of its ``self.config``.

    Examples
    --------
    ::

        class MyWorkflow:
            @with_scratch_temp_dir
            def run(self): ...

    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        with scratch_temp_dir(self.config.run_config.product_path_group.scratch_path):
            return method(self, *args, **kwargs)

    return wrapper  # type: ignore[return-value]
