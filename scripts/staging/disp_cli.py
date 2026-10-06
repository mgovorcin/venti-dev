# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""DISP-S1 CLI for previewing and downloading displacement stack products.

This module provides tools to query, preview, and download OPERA DISP-S1
displacement products for specific frames within date ranges.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import tyro
from opera_utils import get_frame_geodataframe
from opera_utils.disp import search
from opera_utils.disp._download import run_download
from opera_utils.disp._product import DispProductStack
from utils import parse_date

# Constants
DATE_FORMATS = ("%Y-%m-%d", "%Y%m%d")
DEFAULT_NUM_WORKERS = 2
SINGLE_DATE_BUFFER_DAYS = 1
DEFAULT_OUTPUT_DIR = Path("./disp_s1")

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

# Suppress verbose logging from third-party libraries
logging.getLogger("opera_utils").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("asf_search").setLevel(logging.WARNING)


def _normalize_dates(dates: list) -> list[str]:
    """Normalize date objects to ISO date strings.

    Parameters
    ----------
    dates : list
        List of dates as numpy.datetime64, datetime, or string objects.

    Returns
    -------
    list[str]
        Sorted list of unique ISO date strings (YYYY-MM-DD).

    """
    result = []
    for d in dates:
        if isinstance(d, np.datetime64):
            normalized = np.datetime_as_string(d, unit="D")
        elif isinstance(d, datetime):
            normalized = d.date().isoformat()
        else:
            normalized = str(d)
        result.append(normalized)
    return sorted(set(result))


def _validate_frame_in_database(frame_id: int) -> None:
    """Validate that frame ID exists in North America DISP database.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.

    Raises
    ------
    ValueError
        If frame ID is not found in the North America DISP database.

    """
    disp_frame_db = get_frame_geodataframe()
    opera_disp_db = disp_frame_db[disp_frame_db.is_north_america]

    if str(frame_id) not in opera_disp_db.index.astype(str).values:
        msg = (
            f"Frame ID {frame_id} not found in North America DISP database. "
            "Only North American frames are supported."
        )
        raise ValueError(msg)


def _print_product_summary(
    dps: DispProductStack,
    print_urls: bool = False,
    print_ref_dates: bool = False,
    print_sec_dates: bool = False,
    print_dates: bool = False,
) -> None:
    """Print summary of DISP products found.

    Parameters
    ----------
    dps : DispProductStack
        Stack of DISP products.
    frame_id : int
        Frame identifier for logging.
    print_urls : bool, optional
        Whether to print all product URLs. Default is False.
    print_ref_dates : bool, optional
        Whether to print reference dates. Default is False.
    print_sec_dates : bool, optional
        Whether to print secondary dates. Default is False.
    print_dates : bool, optional
        Whether to print all unique dates. Default is False.

    """
    ref_dates = _normalize_dates(dps.reference_dates)
    sec_dates = _normalize_dates(dps.secondary_dates)
    all_dates = sorted(set(ref_dates + sec_dates))

    logger.info(f"Total files: {len(dps.filenames)}")
    logger.info(f"Reference dates: {len(ref_dates)}")
    logger.info(f"Secondary dates: {len(sec_dates)}")
    logger.info(f"Total unique dates: {len(all_dates)}")

    if print_ref_dates:
        logger.info("\nReference dates:")
        for date in ref_dates:
            print(date)

    if print_sec_dates:
        logger.info("\nSecondary dates:")
        for date in sec_dates:
            print(date)

    if print_dates:
        logger.info("\nAll unique dates:")
        for date in all_dates:
            print(date)

    if print_urls:
        logger.info("\nProduct URLs:")
        for url in dps.filenames:
            print(url)


def preview_frame_products(
    frame_id: int,
    start: datetime | None = None,
    end: datetime | None = None,
    print_urls: bool = False,
    print_ref_dates: bool = False,
    print_sec_dates: bool = False,
    print_dates: bool = False,
) -> None:
    """Preview DISP-S1 products available for a frame.

    Queries the OPERA DISP-S1 archive and displays summary statistics
    about available displacement products within the specified date range.
    By default, only counts are shown; specific dates can be printed
    using the optional flags.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    start : datetime or None, optional
        Start date for query. Default is None (no start limit).
    end : datetime or None, optional
        End date for query. Default is None (no end limit).
    print_urls : bool, optional
        Whether to print all product URLs. Default is False.
    print_ref_dates : bool, optional
        Whether to print reference dates list. Default is False.
    print_sec_dates : bool, optional
        Whether to print secondary dates list. Default is False.
    print_dates : bool, optional
        Whether to print all unique dates list. Default is False.

    Raises
    ------
    ValueError
        If frame ID is not in the North America DISP database or if
        no products are found for the specified frame and date range.

    Examples
    --------
    >>> preview_frame_products(frame_id=8887, start=datetime(2024, 1, 1))

    """
    logger.info(f"Querying DISP-S1 products for frame {frame_id}")

    _validate_frame_in_database(frame_id)

    results = search(frame_id, start_datetime=start, end_datetime=end)
    if not results:
        msg = (
            f"No DISP-S1 products found for frame {frame_id} "
            f"in date range [{start}, {end}]"
        )
        raise ValueError(msg)

    dps = DispProductStack(results)
    _print_product_summary(
        dps,
        print_urls,
        print_ref_dates,
        print_sec_dates,
        print_dates,
    )


def _adjust_single_date_range(
    start: datetime | None,
    end: datetime | None,
) -> tuple[datetime | None, datetime | None]:
    """Expand date range if start and end are identical.

    When querying for a single date, expand the range by adding buffer days
    to ensure products are captured.

    Parameters
    ----------
    start : datetime or None
        Start datetime.
    end : datetime or None
        End datetime.

    Returns
    -------
    tuple[datetime or None, datetime or None]
        Adjusted start and end datetimes.

    """
    if start and end and start == end:
        logger.info(
            "Single date query detected. Expanding range by "
            f"±{SINGLE_DATE_BUFFER_DAYS} day(s)"
        )
        buffer = timedelta(days=SINGLE_DATE_BUFFER_DAYS)
        return start - buffer, end + buffer
    return start, end


def download_frame_products(
    frame_id: int,
    output_dir: Path,
    start: datetime | None = None,
    end: datetime | None = None,
    num_workers: int = DEFAULT_NUM_WORKERS,
) -> None:
    """Download DISP-S1 products for a frame.

    Downloads displacement products from the OPERA DISP-S1 archive for
    the specified frame and date range. Products are filtered based on
    the secondary date of each interferogram.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Directory where products will be saved.
    start : datetime or None, optional
        Start date for query (based on secondary date).
        Default is None (no start limit).
    end : datetime or None, optional
        End date for query (based on secondary date).
        Default is None (no end limit).
    num_workers : int, optional
        Number of parallel download workers. Default is 2.

    Raises
    ------
    ValueError
        If frame ID is not in the database or no products are found.

    Notes
    -----
    Date queries are based on the secondary (later) date of each
    interferometric pair. If start and end dates are identical,
    the range is automatically expanded by ±1 day and num_workers
    is set to 1 to ensure the specific product is captured.

    Examples
    --------
    >>> download_frame_products(
    ...     frame_id=8887,
    ...     output_dir=Path("./data"),
    ...     start=datetime(2024, 1, 1),
    ...     end=datetime(2024, 12, 31)
    ... )

    """
    logger.info(f"Downloading DISP-S1 products for frame {frame_id}")

    start_adjusted, end_adjusted = _adjust_single_date_range(start, end)
    workers = 1 if (start and end and start == end) else num_workers

    run_download(
        frame_id,
        start_datetime=start_adjusted,
        end_datetime=end_adjusted,
        output_dir=output_dir,
        num_workers=workers,
    )

    logger.info(f"Download complete: files saved to {output_dir}")


@dataclass
class Preview:
    """CLI for previewing DISP-S1 products.

    Attributes
    ----------
    frame_id : int
        OPERA frame identifier.
    start : str or None
        Start date (YYYY-MM-DD or YYYYMMDD format).
    end : str or None
        End date (YYYY-MM-DD or YYYYMMDD format).
    print_urls : bool
        Whether to print all product URLs.
    print_ref_dates : bool
        Whether to print list of reference dates.
    print_sec_dates : bool
        Whether to print list of secondary dates.
    print_dates : bool
        Whether to print list of all unique dates.

    """

    frame_id: int
    start: str | None = None
    end: str | None = None
    print_urls: bool = False
    print_ref_dates: bool = False
    print_sec_dates: bool = False
    print_dates: bool = False

    def __call__(self) -> None:
        """Execute preview command."""
        preview_frame_products(
            frame_id=self.frame_id,
            start=parse_date(self.start),
            end=parse_date(self.end),
            print_urls=self.print_urls,
            print_ref_dates=self.print_ref_dates,
            print_sec_dates=self.print_sec_dates,
            print_dates=self.print_dates,
        )


@dataclass
class Download:
    """CLI for downloading DISP-S1 products.

    Attributes
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Directory where products will be saved.
    start : str or None
        Start date (YYYY-MM-DD or YYYYMMDD format).
    end : str or None
        End date (YYYY-MM-DD or YYYYMMDD format).
    num_workers : int
        Number of parallel download workers.

    """

    frame_id: int
    output_dir: Path = DEFAULT_OUTPUT_DIR
    start: str | None = None
    end: str | None = None
    num_workers: int = DEFAULT_NUM_WORKERS

    def __call__(self) -> None:
        """Execute download command."""
        download_frame_products(
            frame_id=self.frame_id,
            output_dir=self.output_dir,
            start=parse_date(self.start),
            end=parse_date(self.end),
            num_workers=self.num_workers,
        )


def main() -> None:
    """Run CLI with preview and download subcommands."""
    tyro.extras.subcommand_cli_from_dict(
        {
            "preview": Preview,
            "download": Download,
        }
    )()


if __name__ == "__main__":
    main()
