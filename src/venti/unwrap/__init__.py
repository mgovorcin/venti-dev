# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Unwrap error correction module.

`unwrap_corrections` is the gamma corrector (segments on the displacement,
shifts regions to a common median). `regions` + `cycles` are the v0.5
water-mask estimator (plan T36): whole-cycle shifts per island or cut-off
strip, GNSS-direction veto, gated behind
`CalibrationOptions.unwrap_error_correction` (off until TS-U1).
"""

from .cycles import (
    Decisions,
    RegionDecision,
    apply_shifts,
    estimate_cycles,
    make_unwrap_hook,
    shift_field,
    unwrap_hook_for_pair,
)
from .regions import downsample_labels, largest_region, segment_regions
from .unwrap_corrections import UnwrapCorrector, correct_region_offset

__all__ = [
    "Decisions",
    "RegionDecision",
    "UnwrapCorrector",
    "apply_shifts",
    "correct_region_offset",
    "downsample_labels",
    "estimate_cycles",
    "largest_region",
    "make_unwrap_hook",
    "segment_regions",
    "shift_field",
    "unwrap_hook_for_pair",
]
