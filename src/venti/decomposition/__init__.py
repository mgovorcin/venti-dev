# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""LOS decomposition and projection (PRD section 3.2, R-E3; plan T53).

The `venti[decomposition]` tier: numpy only. `decompose_wls` solves
``[E, U]`` per pixel from two or more look geometries with north fixed from
GNSS; `project_vertical` handles a single geometry with east and north from
GNSS; `decompose` dispatches per pixel and flags the mode.
"""

from .wls import (
    MODE_NONE,
    MODE_PROJECTION,
    MODE_WLS,
    DecompositionResult,
    decompose,
    decompose_wls,
    project_vertical,
)

__all__ = [
    "MODE_NONE",
    "MODE_PROJECTION",
    "MODE_WLS",
    "DecompositionResult",
    "decompose",
    "decompose_wls",
    "project_vertical",
]
