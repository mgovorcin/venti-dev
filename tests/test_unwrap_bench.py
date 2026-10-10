# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""The 14-case F08882 unwrap bench (plan T36.4; trade study ``unwrap_prototype``).

Port of ``harness.py`` + ``score.py`` + ``truth.csv``: on the cached per-pair
inputs (``<bench>/<pair>/cache/{R,G,coh,ws,res}.npy`` + ``meta.yaml``) the
estimator must reproduce the 14 labelled decisions, with no shift on a
region that is not labelled. A regression floor, not proof (2 real errors
among the 14 cases). Set ``VENTI_UNWRAP_BENCH_DIR`` to the ``runs/pairs``
directory to run it; each pair is a full 30 m frame (about 1.2 GB).
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml

from venti.unwrap.cycles import Decisions, estimate_cycles
from venti.workflow.config import UnwrapOptions

BENCH = os.environ.get("VENTI_UNWRAP_BENCH_DIR")
TRUTH = Path(__file__).parent / "data" / "unwrap_truth.csv"
PIXEL_M = 30.0

pytestmark = pytest.mark.skipif(not BENCH, reason="VENTI_UNWRAP_BENCH_DIR not set")


def load_truth() -> list[dict[str, str]]:
    with TRUTH.open() as fh:
        return list(csv.DictReader(fh))


def score(decisions: list[Decisions], truth: list[dict[str, str]]) -> dict[str, Any]:
    """Port of ``score.py``: correct / missed / false_shift / wrong_n / unverified."""
    got = {(d.pair, r.region): r for d in decisions for r in d.rows}
    counts: dict[str, int] = {
        "correct": 0,
        "missed": 0,
        "false_shift": 0,
        "wrong_n": 0,
        "unverified": 0,
    }
    lines: list[str] = []
    for t in truth:
        r = got.get((t["pair"], int(t["region"])))
        cyc = r.cycles if r else 0
        reason = r.reason if r else "not measured"
        if t["expect"] in ("noshift", "0"):
            kind = "correct" if cyc == 0 else "false_shift"
        else:
            want = int(t["expect"])
            kind = "correct" if cyc == want else ("missed" if cyc == 0 else "wrong_n")
        counts[kind] += 1
        jump = f"{r.jump:+.2f}" if r and np.isfinite(r.jump) else "nan"
        lines.append(
            f"{'ok ' if kind == 'correct' else 'BAD'} {t['pair']} r{t['region']:>3} "
            f"expect {t['expect']:>7} jump {jump:>6} -> {reason}"
        )
    labelled = {(t["pair"], int(t["region"])) for t in truth}
    for (pair, region), r in got.items():
        if r.cycles and (pair, region) not in labelled:
            counts["unverified"] += 1
            lines.append(f"??? {pair} r{region} unlabelled -> {r.reason}")
    return {**counts, "lines": lines}


def pair_dirs() -> list[Path]:
    root = Path(BENCH or ".")
    return sorted(p for p in root.iterdir() if (p / "cache" / "ws.npy").exists())


def run_bench(options: UnwrapOptions) -> list[Decisions]:
    decisions = []
    for p in pair_dirs():
        cache = p / "cache"
        meta = yaml.safe_load((cache / "meta.yaml").read_text())
        arrays = {
            n: np.load(cache / f"{n}.npy") for n in ("R", "G", "coh", "ws", "res")
        }
        decisions.append(
            estimate_cycles(
                arrays["R"],
                arrays["ws"],
                arrays["G"],
                arrays["coh"],
                float(meta["half"]),
                PIXEL_M,
                options,
                inversion_residual=arrays["res"],
                pair=str(meta["pair"]),
            )
        )
        del arrays
    return decisions


def test_truth_table_is_the_prototype_one():
    truth = load_truth()
    assert len(truth) == 14
    assert {t["expect"] for t in truth} == {"-3", "-1", "0", "noshift"}


@pytest.mark.parametrize(
    "options",
    [
        pytest.param(UnwrapOptions(), id="ref_gnss_veto"),
        pytest.param(
            UnwrapOptions(residual_gate_cycles=0.05), id="ref_gnss_veto_resid"
        ),
    ],
)
def test_estimator_scores_14_of_14(options):
    truth = load_truth()
    dirs = pair_dirs()
    if len(dirs) < 7:
        pytest.skip(f"bench needs the 7 cached pairs, found {len(dirs)}")
    s = score(run_bench(options), truth)
    report = "\n".join(s["lines"])
    assert s["correct"] == 14, report
    assert s["missed"] == s["false_shift"] == s["wrong_n"] == 0, report
    assert s["unverified"] == 0, report
