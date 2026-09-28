"""Workflow module for InSAR displacement processing.

This module provides workflows for InSAR displacement products:

1. Calibration: Using GNSS reference data to calibrate displacement products
2. Decomposition: Converting LOS displacement to ENU components

APIs:
- Configuration management (YAML-based)
- Object-oriented API: CalibrationWorkflow, DecompositionWorkflow classes
- Functional API: run_workflow, calibrate_timeseries, decompose_timeseries functions
"""

from __future__ import annotations

__all__ = [
    "CalibrationState",
    # Calibration API
    "CalibrationWorkflow",
    "DecompositionState",
    # Decomposition API
    "DecompositionWorkflow",
    # Configuration
    "WorkflowConfig",
    # Workflow types
    "WorkflowType",
    "calibrate_timeseries",
    "create_config_templates",
    "decompose_timeseries",
    # Calibration step functions
    "find_event_mask_file",
    "find_reference_point",
    "load_config",
    "load_los_and_mask",
    "process_displacement_file",
    "run_data_staging",
    "run_data_staging_window",
    # Functional API
    "run_workflow",
    "setup_gnss_reference",
]


# Lazy imports to avoid heavy dependencies on import
def __getattr__(name: str):
    """Lazily import workflow classes and functions on first access."""
    # Configuration
    if name == "WorkflowConfig":
        from .config import WorkflowConfig

        globals()["WorkflowConfig"] = WorkflowConfig
        return WorkflowConfig

    if name == "load_config":
        from .config import load_config

        globals()["load_config"] = load_config
        return load_config

    if name == "create_config_templates":
        from .config import create_config_templates

        globals()["create_config_templates"] = create_config_templates
        return create_config_templates

    # Calibration API: the CalibrationWorkflow class plus the standalone step
    # functions it's built on (setup_gnss_reference, load_los_and_mask, etc.)
    if name in [
        "CalibrationWorkflow",
        "CalibrationState",
        "setup_gnss_reference",
        "load_los_and_mask",
        "find_reference_point",
        "find_event_mask_file",
        "process_displacement_file",
    ]:
        from . import calibration

        for attr in (
            "CalibrationWorkflow",
            "CalibrationState",
            "setup_gnss_reference",
            "load_los_and_mask",
            "find_reference_point",
            "find_event_mask_file",
            "process_displacement_file",
        ):
            globals()[attr] = getattr(calibration, attr)

        return globals()[name]

    # Decomposition API
    if name in [
        "DecompositionWorkflow",
        "DecompositionState",
    ]:
        from .decomposition import (
            DecompositionState,
            DecompositionWorkflow,
        )

        globals()["DecompositionWorkflow"] = DecompositionWorkflow
        globals()["DecompositionState"] = DecompositionState

        return globals()[name]

    # Functional API
    if name == "WorkflowType":
        from .run import WorkflowType

        globals()["WorkflowType"] = WorkflowType
        return WorkflowType

    if name == "run_workflow":
        from .run import run_workflow

        globals()["run_workflow"] = run_workflow
        return run_workflow

    if name == "calibrate_timeseries":
        from .run import calibrate_timeseries

        globals()["calibrate_timeseries"] = calibrate_timeseries
        return calibrate_timeseries

    if name == "decompose_timeseries":
        from .run import decompose_timeseries

        globals()["decompose_timeseries"] = decompose_timeseries
        return decompose_timeseries

    if name == "run_data_staging":
        from .run import run_data_staging

        globals()["run_data_staging"] = run_data_staging
        return run_data_staging

    if name == "run_data_staging_window":
        from .run import run_data_staging_window

        globals()["run_data_staging_window"] = run_data_staging_window
        return run_data_staging_window

    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
