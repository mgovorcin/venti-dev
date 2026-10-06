# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""The v0.5 calibration engine (PRD §3.3, plan T29-T36).

Built next to the gamma `venti.surface` entry point, behind the schema-v2
flags in `calibration_options`, so the gamma algorithm stays reproducible:

- `gaps`: fill masked/water cells before the fit so the surface is continuous
  everywhere and never 0 on masked cells (R-S3);
- `loclin` (T30): local-linear kernel with a physical half-response cutoff;
- `weights` (T31): robust x coherence^p weights;
- `remove_restore` (T32): defo/event areas;
- `two_pass` (T33): robust tie -> unwrap hook -> final surface, with the
  component bookkeeping `calibration == sum(cal_*)`.
"""

from .gaps import base_weights, fill_gaps

__all__ = ["base_weights", "fill_gaps"]
