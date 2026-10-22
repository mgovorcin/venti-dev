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


def test_bundled_table_covers_the_benchmark_and_vlm_frames():
    table = load_frame_table()
    assert DEFAULT_TABLE.exists()
    assert table.version
    # T38 frames, T40 benchmark frames 5-8, the Houston and LA descending VLM partners
    assert table.frames() == [
        7081,
        8622,
        8882,
        8886,
        16940,
        16941,
        18904,
        23211,
        35991,
        36542,
        38238,
    ]
    # PRD D3: Hawaii on the Pacific plate, Puerto Rico on the Caribbean
    plates = {
        f: table.apply(AlgorithmParameters(), f).calibration_options.frame.plate
        for f in (23211, 35991, 36542)
    }
    assert plates == {23211: "PA", 35991: "CA", 36542: "NA"}
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
    # untouched settings keep their defaults; the input is not mutated
    assert cal.surface.regions == "separate"  # T62 edge pass (0.4-draft)
    assert cal.surface.model_copy(update={"regions": "none"}) == (
        CalibrationOptions().surface
    )
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


def test_ts_t1_tropo_modes_in_the_0_3_to_1_5_km_band():
    """TS-T1 (2026-10-14): only humid Puerto Rico gets the correction in the
    band; New York stays off and Montana keeps the default (off)."""
    table = load_frame_table()

    def mode(frame):
        return table.apply(AlgorithmParameters(), frame).calibration_options.tropo.mode

    assert mode(35991) == "stratified"  # Puerto Rico, 0.70 km
    assert mode(8622) == "off"  # New York, 0.42 km
    assert mode(7081) == AlgorithmParameters().calibration_options.tropo.mode
    assert mode(16940) == "stratified"  # Los Angeles, 1.56 km (earlier study)


def test_la_descending_partner_mirrors_la():
    """0.7-draft (2026-10-20): F18904, the LA descending VLM partner, carries
    F16940's island regions and stratified tropo (T57 LA release re-run)."""
    table = load_frame_table()
    asc = table.apply(AlgorithmParameters(), 16940).calibration_options
    desc = table.apply(AlgorithmParameters(), 18904).calibration_options
    assert desc.surface.regions == asc.surface.regions == "separate"
    assert desc.tropo.mode == asc.tropo.mode == "stratified"
    unlisted = table.apply(AlgorithmParameters(), 99999).calibration_options
    assert desc.uncertainty.k_grid == unlisted.uncertainty.k_grid  # no TS-G1 fit


def test_hawaii_takes_stratified_tropo():
    """0.8-draft (2026-10-21): Hawaii (1.62 km) gets stratified tropo; the bias
    rise that failed it is not significant (disp2vlm_validation decision 0004)."""
    table = load_frame_table()
    opts = table.apply(AlgorithmParameters(), 23211).calibration_options
    assert opts.tropo.mode == "stratified"
    assert opts.frame.plate == "PA"
