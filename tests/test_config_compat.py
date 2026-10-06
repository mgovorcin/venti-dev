# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Compatibility tests between cal-disp's algorithm_parameters.yaml and Venti.

cal-disp (the operational DISP-CAL SAS) ships `algorithm_parameters.yaml`
files whose `calibration_options` block overlaps Venti's `CalibrationOptions`.
The gamma 0.3 golden configuration is vendored under tests/data so a Venti
change that breaks loading it, or silently starts ignoring a key, fails here
before cal-disp notices.

Known gap (plan T19): Venti's `CalibrationOptions` does not forbid extra keys,
and cal-disp's downsampling keys are passed to `estimate_calibration_surface`
as function arguments rather than modelled. The second test pins that set so
it can only change on purpose.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from venti.workflow.config import AlgorithmParameters, CalibrationOptions

GAMMA_YAML = Path(__file__).parent / "data/caldisp_algorithm_parameters_gamma.yaml"

# cal-disp keys that Venti's options model does not carry (they are arguments
# of estimate_calibration_surface in cal-disp's call). Modelling them is T19.
KEYS_NOT_MODELLED_BY_VENTI = {
    "downsample_factor",
    "downsample_method",
    "downsample_weighted",
}


@pytest.fixture(scope="module")
def gamma_yaml() -> dict:
    with GAMMA_YAML.open() as f:
        return yaml.safe_load(f)


def test_gamma_calibration_options_load(gamma_yaml):
    opts = CalibrationOptions(**gamma_yaml["calibration_options"])
    # The values that define the gamma 0.3 algorithm must survive the round trip.
    assert opts.grid_type == "constant"
    assert opts.reference_frame == "IGS20"
    assert opts.unwrap_error_correction is False
    assert opts.apply_tropo_correction is True
    assert opts.apply_solid_earth_tide_correction is True
    assert opts.window_size_meters == pytest.approx(600000.0)
    assert opts.posting_meters == pytest.approx(30.0)
    assert opts.mask_fit_residual_outliers is True
    assert opts.weight_fit_by_gnss_uncertainty is False
    assert opts.calibration_surface_smoothing_sigma == pytest.approx(0.0)
    assert opts.residual_outlier_mad_threshold is None


def test_gamma_keys_venti_does_not_model_are_exactly_the_known_ones(gamma_yaml):
    gamma_keys = set(gamma_yaml["calibration_options"])
    modelled = set(CalibrationOptions.model_fields)
    assert gamma_keys - modelled == KEYS_NOT_MODELLED_BY_VENTI
    # And every key Venti does model that gamma sets is accepted (no typos).
    assert gamma_keys & modelled == gamma_keys - KEYS_NOT_MODELLED_BY_VENTI


def test_gamma_file_loads_as_algorithm_parameters(tmp_path):
    params = AlgorithmParameters.from_yaml(GAMMA_YAML)
    assert params.calibration_options.window_size_meters == pytest.approx(600000.0)
    # Round trip: write and re-read gives the same calibration options.
    out = tmp_path / "roundtrip.yaml"
    params.to_yaml(out)
    again = AlgorithmParameters.from_yaml(out)
    assert again.calibration_options == params.calibration_options
