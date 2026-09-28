"""Main CLI entry point for Venti.

This module provides the command-line interface for the Venti calibration workflow.

Usage:
    venti config --output-dir DIR                                # Config templates
    venti run --config-file runconfig.yaml                       # All products
    venti run-single --config-file runconfig.yaml --disp-file epoch.nc  # One

Log messages go to the console and, if the runconfig sets ``log_file``, to
that file as well (see `venti.log_setup.configure_logging`).
"""

from __future__ import annotations

import logging
import sys
from enum import StrEnum
from typing import Literal

import tyro

from .log_setup import configure_logging

# Not __name__: under `python -m venti` that is "__main__", outside the
# "venti" logger that configure_logging sets up.
logger = logging.getLogger("venti.cli")

LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]


class Command(StrEnum):
    """Available Venti commands."""

    config = "config"
    run = "run"
    run_single = "run-single"


def config_command(output_dir: str = ".") -> None:
    """Generate configuration file templates.

    Generates two files:
    - runconfig.yaml: Run-specific settings with BOTH calibration and decomposition
      input groups (fill in the one you need based on workflow type)
    - algorithm_parameters.yaml: Algorithm defaults (GNSS, decomposition settings)

    The runconfig.yaml is structured similar to algorithm_parameters.yaml:
    - calibration_input_group section for calibration workflow
    - decomposition_input_group section for decomposition workflow
    - Users only fill in the section they need based on workflow_name

    Parameters
    ----------
    output_dir : str
        Directory to write configuration files (default: current directory)

    Examples
    --------
    ::

        # Generate config files in current directory
        python -m venti config

        # Generate config files in specific directory
        python -m venti config --output-dir configs/

    """
    from .workflow.config import create_config_templates

    configure_logging()
    try:
        runconfig_path, params_path = create_config_templates(output_dir)
        logger.info("Configuration templates created:")
        logger.info(f"  1. {runconfig_path}")
        logger.info(f"  2. {params_path}")
        logger.info("")
        logger.info("The runconfig.yaml contains TWO input group sections:")
        logger.info("  - calibration_input_group (for calibration workflow)")
        logger.info("  - decomposition_input_group (for decomposition workflow)")
        logger.info("")
        logger.info("Next steps:")
        logger.info("  1. Edit runconfig.yaml:")
        logger.info("     - Fill in calibration_input_group if doing calibration")
        logger.info("     - Fill in decomposition_input_group if doing decomposition")
        logger.info("     - Set workflow_name to 'calibrate' or 'decompose'")
        logger.info("  2. Review algorithm_parameters.yaml (defaults usually work)")
        logger.info("  3. Run: venti run --config-file runconfig.yaml")
    except Exception:
        logger.exception("Failed to create configuration templates")
        sys.exit(1)


def run_command(
    config_file: str, n_workers: int = 1, log_level: LogLevel = "INFO"
) -> None:
    """Run the Venti calibration workflow.

    Parameters
    ----------
    config_file : str
        Path to YAML configuration file.
    n_workers : int, optional
        Displacement files processed concurrently in threads, by default 1.
        I/O overlaps across files; the surface fits run one at a time, each
        using all CPUs. 2-4 is a good range.
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR), default: INFO.

    Examples
    --------
    ::

        venti run --config-file runconfig.yaml
        venti run --config-file runconfig.yaml --n-workers 4
        venti run --config-file runconfig.yaml --log-level DEBUG

    """
    from .workflow.calibration import CalibrationWorkflow
    from .workflow.config import load_config

    configure_logging(level=log_level)

    try:
        # Load configuration
        logger.info(f"Loading configuration from: {config_file}")
        config = load_config(config_file)
        configure_logging(log_file=config.run_config.log_file)
        logger.info("Configuration loaded successfully")
        logger.info(f"  Input directory: {config.input_options.input_files}")
        logger.info(
            f"  Output directory: {config.run_config.product_path_group.product_path}"
        )
        logger.info(f"  Grid type: {config.grid_settings.grid_type}")
        logger.info(f"  Reference frame: {config.grid_settings.reference_frame}")

        # Run workflow using OO API
        logger.info("Starting calibration workflow...")
        workflow = CalibrationWorkflow(config=config)
        workflow.run(n_workers=n_workers)

        logger.info("Workflow completed successfully!")

    except FileNotFoundError:
        logger.exception("File not found")
        sys.exit(1)
    except ValueError:
        logger.exception("Configuration error")
        sys.exit(1)
    except Exception:
        logger.exception("Workflow failed")
        sys.exit(1)


def run_single_command(
    config_file: str,
    disp_file: str,
    tropo_ref_file: str | None = None,
    tropo_sec_file: str | None = None,
    log_level: LogLevel = "INFO",
) -> None:
    """Calibrate a single displacement file.

    Parameters
    ----------
    config_file : str
        Path to YAML configuration file.
    disp_file : str
        Path to the NetCDF displacement file to calibrate.
    tropo_ref_file : str, optional
        Per-epoch tropospheric correction GeoTIFF for the reference date.
    tropo_sec_file : str, optional
        Per-epoch tropospheric correction GeoTIFF for the secondary date.
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR), default: INFO.

    Examples
    --------
    ::

        venti run-single --config-file runconfig.yaml --disp-file epoch_001.nc
        venti run-single --config-file runconfig.yaml --disp-file epoch_001.nc \
            --tropo-ref-file /data/tropo/ref.tif --tropo-sec-file /data/tropo/sec.tif
        venti run-single --config-file runconfig.yaml --disp-file epoch_001.nc \
            --log-level DEBUG

    """
    from pathlib import Path

    from .workflow.calibration import CalibrationWorkflow
    from .workflow.config import load_config

    configure_logging(level=log_level)

    try:
        logger.info(f"Loading configuration from: {config_file}")
        config = load_config(config_file)
        configure_logging(log_file=config.run_config.log_file)
        logger.info("Starting single-file calibration...")
        workflow = CalibrationWorkflow(config=config)
        state = workflow.run_single(
            disp_file=Path(disp_file),
            tropo_ref_file=Path(tropo_ref_file) if tropo_ref_file is not None else None,
            tropo_sec_file=Path(tropo_sec_file) if tropo_sec_file is not None else None,
        )

        if state.n_files_failed:
            logger.error("Calibration failed — check logs above.")
            sys.exit(1)

        logger.info(f"Output: {state.output_files[0]}")

    except FileNotFoundError:
        logger.exception("File not found")
        sys.exit(1)
    except ValueError:
        logger.exception("Configuration error")
        sys.exit(1)
    except Exception:
        logger.exception("Workflow failed")
        sys.exit(1)


def main() -> None:
    """Run the main CLI with subcommands."""
    # Use tyro for CLI parsing with subcommands
    tyro.extras.subcommand_cli_from_dict(
        {
            Command.config: config_command,
            Command.run: run_command,
            Command.run_single: run_single_command,
        }
    )


if __name__ == "__main__":
    main()
