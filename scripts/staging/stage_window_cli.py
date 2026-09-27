r"""CLI for staging all ancillary data for a frame over a time window.

Downloads every DISP-S1 product whose secondary date falls within the given
window and stages all ancillary data required for downstream calibration.
Frame-level assets (DEM, LOS geometry, GNSS velocities) are produced once and
shared across all products.  Tropospheric corrections are batched over all
unique epoch sensing times and combined into per-product differential files.

Output directory structure::

    <output_dir>/
    ├── disp_s1/   # All downloaded DISP-S1 NetCDF files
    ├── dem/        # Shared GLO-30 DEM
    ├── los/        # Shared LOS ENU and incidence angle rasters
    ├── tropo/      # Per-product differential tropo corrections
    └── gnss/       # Shared UNR GNSS velocities

Examples
--------
Stage all products for frame 8887 in 2016::

    python stage_window_cli.py --frame-id 8887 --start 2016-01-01 --end 2016-12-31

Stage without tropospheric corrections::

    python stage_window_cli.py --frame-id 8887 --start 2016-01-01 --end 2016-12-31 \
        --skip-tropo

Custom output directory and parallelism::

    python stage_window_cli.py --frame-id 8887 --start 2016-01-01 --end 2016-12-31 \
        --output-dir /data/opera/frame_8887 --num-workers 8

"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import tyro

from venti.workflow import run_data_staging_window


@dataclass
class StageWindowCLI:
    """Stage all ancillary data for a frame over a secondary-date time window.

    Attributes
    ----------
    frame_id : int
        OPERA frame identifier.
    start : str
        Start of the secondary-date window (YYYY-MM-DD or YYYYMMDD).
    end : str
        End of the secondary-date window (YYYY-MM-DD or YYYYMMDD).
    output_dir : Path
        Root directory for all staged outputs.
    num_workers : int
        Number of parallel workers for downloads and processing.
    dem_buffer : float
        Buffer in metres around the frame extent for DEM generation.
    skip_tropo : bool
        Skip tropospheric correction processing.
    skip_gnss : bool
        Skip UNR GNSS download and velocity estimation.
    gnss_reference_frame : str
        GNSS reference frame for UNR data ('IGS20' or 'IGS14').
    gnss_padding : float
        Extra padding in metres beyond frame bounds when searching for
        GNSS stations.
    gnss_start_year : float
        Exclude GNSS observations before this decimal year when estimating
        velocities.
    log_file : Path or None
        Also write log messages to this file.

    """

    frame_id: int
    start: str
    end: str
    output_dir: Path = Path("./staging")
    num_workers: int = 4
    dem_buffer: float = 10_000.0
    skip_tropo: bool = False
    skip_gnss: bool = False
    gnss_reference_frame: str = "IGS20"
    gnss_padding: float = 0.0
    gnss_start_year: float = 2014.0
    log_file: Path | None = None

    def __call__(self) -> None:
        """Execute the window staging workflow."""
        run_data_staging_window(
            frame_id=self.frame_id,
            start=self.start,
            end=self.end,
            output_dir=self.output_dir,
            num_workers=self.num_workers,
            dem_buffer=self.dem_buffer,
            skip_tropo=self.skip_tropo,
            skip_gnss=self.skip_gnss,
            gnss_reference_frame=self.gnss_reference_frame,
            gnss_padding=self.gnss_padding,
            gnss_start_year=self.gnss_start_year,
            log_file=self.log_file,
        )


def main() -> None:
    """Run CLI."""
    tyro.cli(StageWindowCLI)()


if __name__ == "__main__":
    main()
