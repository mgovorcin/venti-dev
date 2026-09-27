r"""CLI for staging all ancillary data for a single OPERA DISP-S1 frame.

Runs the full data staging workflow:
  1. Download DISP-S1 product
  2. Generate GLO30 DEM
  3. Generate LOS ENU and incidence angle rasters
  4. Download and apply tropospheric corrections
  5. Download UNR GNSS station timeseries and compute velocities

Examples
--------
Stage all data for frame 8887 on 2016-06-15::

    python stage_frame_cli.py --frame-id 8887 --date 2016-06-15

Stage without tropospheric corrections or GNSS::

    python stage_frame_cli.py --frame-id 8887 --date 2016-06-15 \\
        --skip-tropo --skip-gnss

Custom output directory and parallelism::

    python stage_frame_cli.py --frame-id 8887 --date 2016-06-15 \\
        --output-dir /data/opera/frame_8887 --num-workers 8

"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import tyro

from venti.workflow import run_data_staging


@dataclass
class StageFrameCLI:
    """Stage all ancillary data for a single OPERA DISP-S1 frame.

    Attributes
    ----------
    frame_id : int
        OPERA frame identifier.
    date : str
        Secondary date of the interferogram to stage (YYYY-MM-DD or YYYYMMDD).
    output_dir : Path
        Root directory for all staged outputs.
    num_workers : int
        Number of parallel workers for downloads and processing.
    dem_buffer : float
        Buffer in meters around the frame extent for DEM generation.
    skip_tropo : bool
        Skip tropospheric correction processing.
    skip_gnss : bool
        Skip UNR GNSS download and velocity estimation.
    gnss_reference_frame : str
        GNSS reference frame for UNR data ('IGS20' or 'IGS14').
    gnss_padding : float
        Extra padding in meters beyond frame bounds when searching for
        GNSS stations.
    gnss_start_year : float
        Exclude GNSS observations before this decimal year when estimating
        velocities.
    log_file : Path or None
        Also write log messages to this file.

    """

    frame_id: int
    date: str
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
        """Execute the staging workflow."""
        run_data_staging(
            frame_id=self.frame_id,
            date=self.date,
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
    tyro.cli(StageFrameCLI)()


if __name__ == "__main__":
    main()
