# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Frame-parameter table: per-frame algorithm overrides (plan T27, PRD §4.3).

The table says, per OPERA frame, which plate the ``plate_motion`` layer uses,
which tropospheric mode applies, the GNSS-sigma inflation ``k`` from trade
study TS-G1, and the CalVal benchmark category. It is a versioned JSON file
with the same shape as cal-disp's ``algorithm_parameters_overrides_json``::

    {
      "version": "1.0",
      "default": {"calibration_options": {"frame": {"plate": "NA"}}},
      "data": {
        "08882": {"calibration_options": {"tropo": {"mode": "off"}, ...}},
        ...
      }
    }

Precedence, lowest first: the algorithm-parameters defaults, the table's
``default`` block, the frame's entry, and finally explicit overrides a caller
passes. cal-disp reads only ``data[frame_id]`` (its runconfig is frozen), so
`FrameParameterTable.materialized` folds ``default`` into every frame entry
to produce a file it can consume as-is.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from .workflow.config import AlgorithmParameters

__all__ = ["DEFAULT_TABLE", "FrameParameterTable", "deep_merge", "load_frame_table"]

DEFAULT_TABLE = Path(__file__).parent / "data" / "frame_parameters.json"


def deep_merge(base: Mapping[str, Any], new: Mapping[str, Any]) -> dict[str, Any]:
    """Return `base` updated by `new`, merging nested mappings (no mutation)."""
    out: dict[str, Any] = deepcopy(dict(base))
    for key, value in new.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = deepcopy(value)
    return out


def _frame_key(frame_id: int | str) -> str:
    """Frames are keyed as zero-padded five-digit strings ("08882")."""
    return f"{int(frame_id):05d}"


class FrameParameterTable(BaseModel):
    """Versioned per-frame algorithm overrides."""

    version: str = Field(description="Table version, recorded in product metadata")
    description: str = Field(
        "", description="What this table encodes and where values come from"
    )
    default: dict[str, Any] = Field(
        default_factory=dict,
        description="Overrides applied to every frame before the frame's own entry",
    )
    data: dict[str, dict[str, Any]] = Field(
        default_factory=dict,
        description="Frame id (zero-padded, e.g. '08882') -> overrides for that frame",
    )

    model_config = {"extra": "forbid"}

    # -- lookup ------------------------------------------------------------------

    def frames(self) -> list[int]:
        """Frame ids with an explicit entry."""
        return sorted(int(k) for k in self.data)

    def has_frame(self, frame_id: int | str) -> bool:
        """Return True if the frame has its own entry (not just the defaults)."""
        return _frame_key(frame_id) in self.data

    def overrides_for(
        self, frame_id: int | str, extra: Mapping[str, Any] | None = None
    ) -> dict[str, Any]:
        """Return the merged overrides for a frame: default < frame entry < `extra`."""
        merged = deep_merge(self.default, self.data.get(_frame_key(frame_id), {}))
        if extra:
            merged = deep_merge(merged, extra)
        return merged

    def apply(
        self,
        params: AlgorithmParameters,
        frame_id: int | str,
        extra: Mapping[str, Any] | None = None,
    ) -> AlgorithmParameters:
        """Return `params` with this frame's overrides applied and validated.

        Override keys are option groups (``calibration_options``, ...); as in
        cal-disp, bare keys are a shorthand for ``calibration_options``.
        """
        overrides = self.overrides_for(frame_id, extra)
        if not overrides:
            return params
        groups = set(type(params).model_fields)
        shaped: dict[str, Any] = {}
        for key, value in overrides.items():
            if key in groups:
                shaped[key] = value
            else:
                shaped.setdefault("calibration_options", {})[key] = value
        data = deep_merge(params.model_dump(), shaped)
        return type(params)(**data)

    # -- export ------------------------------------------------------------------

    def materialized(self) -> dict[str, Any]:
        """Return the table as cal-disp reads it: ``default`` folded into every entry.

        cal-disp's loader returns ``data[str(frame_id)]`` and ignores the other
        keys, so the defaults must be baked into each frame for it to see them.
        """
        return {
            "version": self.version,
            "description": self.description,
            "data": {
                key: deep_merge(self.default, entry) for key, entry in self.data.items()
            },
        }

    def write(self, path: Path | str, *, materialize: bool = False) -> Path:
        """Write the table (optionally materialized for cal-disp) as JSON."""
        path = Path(path)
        payload = self.materialized() if materialize else self.model_dump()
        path.write_text(json.dumps(payload, indent=2) + "\n")
        return path


def load_frame_table(path: Path | str | None = None) -> FrameParameterTable:
    """Load a frame-parameter table; the bundled one when `path` is None."""
    with Path(path or DEFAULT_TABLE).open() as f:
        return FrameParameterTable(**json.load(f))
