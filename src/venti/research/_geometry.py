# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Geometry helpers shared by the research database builders."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

__all__ = ["reproject"]


def reproject(func: Callable[[Any, Any], tuple[Any, Any]], geom: Any) -> Any:
    """Apply a pyproj ``Transformer.transform`` to a shapely geometry.

    Replaces the deprecated ``shapely.ops.transform(func, geom)``:
    `shapely.transform` hands the coordinates over as one (N, 2) array, while
    the pyproj function takes x and y separately.
    """
    from shapely import transform

    return transform(geom, lambda xy: np.column_stack(func(xy[:, 0], xy[:, 1])))
