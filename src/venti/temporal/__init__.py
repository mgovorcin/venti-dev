# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Temporal resampling of displacement time series (PRD D15; plan T54).

Moving-window robust linear fits evaluated at target epochs, so ascending
and descending stacks can be sampled on common dates before decomposition.
"""

from .resample import (
    ResampleMethod,
    resample_dataarray,
    resample_timeseries,
)

__all__ = ["ResampleMethod", "resample_dataarray", "resample_timeseries"]
