# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""algorithm_parameters schema v2 (plan T19): additive, versioned, gamma-equivalent defaults."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from venti.workflow.config import (
    ALGORITHM_SCHEMA_VERSION,
    AlgorithmParameters,
    CalibrationOptions,
)

GAMMA_YAML = Path(__file__).parent / "data/caldisp_algorithm_parameters_gamma.yaml"


def test_defaults_reproduce_the_gamma_algorithm():
    """Every new option defaults to the legacy behaviour (nothing changes on upgrade)."""
    opts = CalibrationOptions()
    assert opts.surface.method == "windowed_plane"
    assert opts.surface.fill_gaps is False
    assert opts.surface.two_pass is False
    assert opts.weights.coherence_power == 0.0
    assert opts.weights.robust is False
    assert opts.tropo.mode == "legacy"
    assert opts.gnss.buffer_meters == 0.0
    assert opts.gnss.exclude_defo_nodes is False
    assert opts.gnss.reprocessing is False
    assert opts.uncertainty.k_grid == 1.0
    assert opts.unwrap_error_correction is False  # R-U1
    assert opts.unwrap.whole_cycles_only is True
    assert opts.unwrap.free_offsets is False
    assert opts.downsample_factor == 1


def test_gamma_file_is_version_1_and_upgrades_in_memory(caplog):
    with GAMMA_YAML.open() as f:
        raw = yaml.safe_load(f)
    assert "schema_version" not in raw
    with caplog.at_level("INFO"):
        params = AlgorithmParameters.from_yaml(GAMMA_YAML)
    assert params.schema_version == ALGORITHM_SCHEMA_VERSION == 2
    assert "upgrading" in caplog.text
    cal = params.calibration_options
    # the gamma values survive, including the keys Venti did not model before
    assert cal.window_size_meters == pytest.approx(600000.0)
    assert cal.downsample_factor == 6
    assert cal.downsample_method == "mean"
    assert cal.downsample_weighted is False
    assert cal.unwrap_error_correction is False
    # and the new groups are at their gamma-equivalent defaults
    assert cal.surface == CalibrationOptions().surface
    assert cal.tropo.mode == "legacy"


def test_round_trip_keeps_version_and_values(tmp_path):
    params = AlgorithmParameters.from_yaml(GAMMA_YAML)
    params.calibration_options.surface.method = "loclin"
    params.calibration_options.weights.coherence_power = 8.0
    out = tmp_path / "v2.yaml"
    params.to_yaml(out)
    with out.open() as f:
        raw = yaml.safe_load(f)
    assert raw["schema_version"] == 2
    assert raw["calibration_options"]["surface"]["method"] == "loclin"
    again = AlgorithmParameters.from_yaml(out)
    assert again == params


def test_unknown_keys_are_rejected():
    with pytest.raises(ValidationError, match="not_a_real_option"):
        CalibrationOptions(not_a_real_option=1)
    with pytest.raises(ValidationError, match="typo"):
        CalibrationOptions(surface={"typo": 1})


def test_unsupported_schema_version():
    with pytest.raises(ValueError, match="schema_version 3"):
        AlgorithmParameters.from_dict({"schema_version": 3})


@pytest.mark.parametrize(
    ("group", "field", "bad"),
    [
        ("surface", "cutoff_wavelength_meters", 0),
        ("surface", "method", "spline"),
        ("weights", "coherence_power", -1),
        ("weights", "filled_pixel_weight", 2),
        ("tropo", "mode", "sometimes"),
        ("gnss", "buffer_meters", -5),
        ("unwrap", "min_region_area", 0),
        ("unwrap", "region_source", "phase"),
    ],
)
def test_ranges_and_literals_are_validated(group, field, bad):
    with pytest.raises(ValidationError):
        CalibrationOptions(**{group: {field: bad}})


def test_k_grid_accepts_number_or_frame_table():
    assert CalibrationOptions(uncertainty={"k_grid": 3.9}).uncertainty.k_grid == 3.9
    assert (
        CalibrationOptions(uncertainty={"k_grid": "frame_table"}).uncertainty.k_grid
        == "frame_table"
    )
    with pytest.raises(ValidationError):
        CalibrationOptions(uncertainty={"k_grid": "per_pixel"})


def test_written_defaults_load_back(tmp_path):
    out = tmp_path / "defaults.yaml"
    AlgorithmParameters().to_yaml(out)
    text = out.read_text()
    assert "schema_version: 2" in text
    assert AlgorithmParameters.from_yaml(out) == AlgorithmParameters()
