# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Frame-parameter table (plan T27)."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from venti.frames import (
    DEFAULT_TABLE,
    FrameParameterTable,
    deep_merge,
    load_frame_table,
)
from venti.workflow.config import AlgorithmParameters, CalibrationOptions


def test_bundled_table_loads_and_covers_the_first_four_frames():
    table = load_frame_table()
    assert DEFAULT_TABLE.exists()
    assert table.version
    assert table.frames() == [8622, 8882, 8886, 16940]
    assert table.has_frame(8882)
    assert table.has_frame("08882")
    assert not table.has_frame(1)


def test_precedence_default_then_frame_then_extra():
    table = FrameParameterTable(
        version="t",
        default={
            "calibration_options": {
                "frame": {"plate": "NA"},
                "uncertainty": {"k_grid": 2.0},
            }
        },
        data={"00001": {"calibration_options": {"frame": {"plate": "PA"}}}},
    )
    merged = table.overrides_for(1)
    assert (
        merged["calibration_options"]["frame"]["plate"] == "PA"
    )  # frame beats default
    assert merged["calibration_options"]["uncertainty"]["k_grid"] == 2.0  # default kept
    merged = table.overrides_for(
        1, extra={"calibration_options": {"uncertainty": {"k_grid": 5.0}}}
    )
    assert (
        merged["calibration_options"]["uncertainty"]["k_grid"] == 5.0
    )  # extra beats all
    # unknown frame -> defaults only
    assert table.overrides_for(99) == table.default


def test_apply_to_algorithm_parameters_from_the_bundled_table():
    table = load_frame_table()
    params = table.apply(AlgorithmParameters(), 16940)
    cal = params.calibration_options
    assert cal.frame.plate == "NA"
    assert cal.frame.name == "Los Angeles"
    assert cal.tropo.mode == "stratified"
    assert cal.uncertainty.k_grid == 2.69  # TS-G1 fit
    # a frame without a fit keeps the default k
    assert (
        table.apply(AlgorithmParameters(), 99).calibration_options.uncertainty.k_grid
        == 3.9
    )
    # untouched groups keep their defaults; the input is not mutated
    assert cal.surface == CalibrationOptions().surface
    assert AlgorithmParameters().calibration_options.frame.plate == "NA"
    houston = table.apply(AlgorithmParameters(), "08882").calibration_options
    assert houston.tropo.mode == "off"
    assert houston.frame.benchmark_category == "coastal_subsidence_islands"


def test_shorthand_keys_target_calibration_options():
    table = FrameParameterTable(
        version="t", data={"00007": {"tropo": {"mode": "full"}}}
    )
    params = table.apply(AlgorithmParameters(), 7)
    assert params.calibration_options.tropo.mode == "full"


def test_invalid_override_is_rejected_on_apply():
    table = FrameParameterTable(
        version="t",
        data={"00007": {"calibration_options": {"tropo": {"mode": "maybe"}}}},
    )
    with pytest.raises(ValidationError):
        table.apply(AlgorithmParameters(), 7)


def test_plate_validation():
    assert CalibrationOptions(frame={"plate": "pa"}).frame.plate == "PA"
    assert CalibrationOptions(frame={"plate": "CARB"}).frame.plate == "CARB"
    with pytest.raises(ValidationError, match="Unknown plate"):
        CalibrationOptions(frame={"plate": "XX"})


def test_materialized_is_what_cal_disp_reads(tmp_path):
    """cal-disp's loader returns data[str(frame_id)] and ignores 'default'."""
    table = load_frame_table()
    out = table.write(tmp_path / "overrides.json", materialize=True)
    raw = json.loads(out.read_text())
    assert "default" not in raw
    entry = raw["data"][
        "08882"
    ]  # exactly what cal-disp's _parse_algorithm_overrides returns
    assert entry["calibration_options"]["frame"]["plate"] == "NA"  # default folded in
    assert entry["calibration_options"]["uncertainty"]["k_grid"] == 1.35  # TS-G1
    assert entry["calibration_options"]["tropo"]["mode"] == "off"
    # and applying that entry alone reproduces the full apply()
    via_entry = FrameParameterTable(version="t", data={"08882": entry}).apply(
        AlgorithmParameters(), 8882
    )
    assert via_entry == table.apply(AlgorithmParameters(), 8882)


def test_write_round_trip(tmp_path):
    table = load_frame_table()
    out = table.write(tmp_path / "table.json")
    assert load_frame_table(out) == table


def test_deep_merge_does_not_mutate():
    base = {"a": {"b": 1, "c": 2}}
    new = {"a": {"b": 9}, "d": 4}
    merged = deep_merge(base, new)
    assert merged == {"a": {"b": 9, "c": 2}, "d": 4}
    assert base == {"a": {"b": 1, "c": 2}}


def test_table_rejects_unknown_top_level_keys():
    with pytest.raises(ValidationError):
        FrameParameterTable(version="t", frames={})
