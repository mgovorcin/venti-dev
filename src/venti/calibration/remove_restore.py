# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Remove-restore: curated defo and event areas (PRD A5, R-S5, §2.8; plan T32).

The calibration field is long-wavelength secular GNSS motion only. Fast local
deformation (subsiding basins, coseismic steps, volcanic episodes) is
*excluded from the fit* and the surface is interpolated through it from the
surrounding ground, so InSAR keeps the full extent of the deforming area and
the GNSS tie is not pulled by it (ADR-0009).

Areas are curated offline as versioned GeoJSON (``defo_area_db_json`` and
``event_db_json`` in the frozen cal-disp runconfig):

- **defo areas**: ``id``, ``name``, polygon; persistent.
- **events**: the same plus ``t0`` (ISO date), optional ``t1`` (end of the
  window during which pairs are affected, e.g. a postseismic phase),
  ``magnitude``, ``source``. An event is applied to a pair only when the pair
  spans ``t0`` (or overlaps ``[t0, t1]`` when ``t1`` is given).

Both files carry a top-level ``version`` (a foreign member of the
FeatureCollection), recorded in the product provenance.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field, field_validator

from ..spatial.interpolation import GridLike, grid_coordinates

logger = logging.getLogger(__name__)

__all__ = [
    "Area",
    "AreaDB",
    "Event",
    "EventDB",
    "RemoveRestore",
    "load_area_db",
    "load_event_db",
    "remove_restore_mask",
    "sigma_inflation_inside",
]


def _as_datetime(value: Any) -> datetime:
    if isinstance(value, datetime):
        return value
    if isinstance(value, date):
        return datetime(value.year, value.month, value.day)
    text = str(value).replace("Z", "+00:00")
    return datetime.fromisoformat(text).replace(tzinfo=None)


class Area(BaseModel):
    """A persistent deforming area (basin, volcano) excluded from the fit."""

    id: str
    name: str = ""
    geometry: dict[str, Any] = Field(description="GeoJSON geometry in EPSG:4326")
    source: str = ""

    model_config = {"extra": "allow"}

    def shape(self) -> Any:
        """Shapely geometry (requires shapely, part of the core)."""
        from shapely.geometry import shape

        return shape(self.geometry)


class Event(Area):
    """A transient event (earthquake, eruption) with a time window."""

    t0: datetime = Field(description="Event time; a pair is affected when it spans t0")
    t1: datetime | None = Field(
        None, description="End of the affected window (postseismic); default none"
    )
    magnitude: float | None = None

    @field_validator("t0", "t1", mode="before")
    @classmethod
    def _parse_time(cls, value: Any) -> datetime | None:
        return None if value is None else _as_datetime(value)

    def affects(self, reference_date: datetime, secondary_date: datetime) -> bool:
        """Return True if the pair spans `t0` (or overlaps ``[t0, t1]``)."""
        ref, sec = sorted((_as_datetime(reference_date), _as_datetime(secondary_date)))
        if self.t1 is None:
            return ref <= self.t0 <= sec
        return ref <= self.t1 and sec >= self.t0


class _DB(BaseModel):
    version: str
    description: str = ""

    model_config = {"extra": "allow"}

    @classmethod
    def _from_geojson(
        cls, data: dict[str, Any], item_cls: type[Area]
    ) -> dict[str, Any]:
        if data.get("type") != "FeatureCollection":
            msg = f"expected a GeoJSON FeatureCollection, got type={data.get('type')!r}"
            raise ValueError(msg)
        version = data.get("version") or data.get("properties", {}).get("version")
        if not version:
            msg = (
                "the GeoJSON needs a top-level 'version' "
                "(recorded in the product provenance)"
            )
            raise ValueError(msg)
        items = []
        for feat in data.get("features", []):
            props = dict(feat.get("properties") or {})
            props["geometry"] = feat["geometry"]
            if "id" not in props and "id" in feat:
                props["id"] = str(feat["id"])
            items.append(item_cls(**props))
        return {
            "version": str(version),
            "description": data.get("description", ""),
            "items": items,
        }


class AreaDB(_DB):
    """Curated persistent deforming areas."""

    items: list[Area] = Field(default_factory=list)

    @classmethod
    def from_geojson(cls, data: dict[str, Any]) -> AreaDB:
        return cls(**cls._from_geojson(data, Area))

    def shapes(self) -> list[Any]:
        return [a.shape() for a in self.items]


class EventDB(_DB):
    """Curated transient events."""

    items: list[Event] = Field(default_factory=list)

    @classmethod
    def from_geojson(cls, data: dict[str, Any]) -> EventDB:
        return cls(**cls._from_geojson(data, Event))

    def active(self, reference_date: datetime, secondary_date: datetime) -> list[Event]:
        """Events that affect the pair."""
        return [e for e in self.items if e.affects(reference_date, secondary_date)]


def load_area_db(path: Path | str) -> AreaDB:
    """Load a defo-area GeoJSON."""
    with Path(path).open() as f:
        return AreaDB.from_geojson(json.load(f))


def load_event_db(path: Path | str) -> EventDB:
    """Load an event GeoJSON."""
    with Path(path).open() as f:
        return EventDB.from_geojson(json.load(f))


# --------------------------------------------------------------------------


def _pixel_centres_lonlat(
    grid: GridLike, crs_epsg: int
) -> tuple[np.ndarray, np.ndarray]:
    from pyproj import Transformer

    gx, gy = grid_coordinates(grid)
    xx, yy = np.meshgrid(np.asarray(gx, float), np.asarray(gy, float))
    if crs_epsg == 4326:
        return xx, yy
    to_ll = Transformer.from_crs(f"EPSG:{crs_epsg}", "EPSG:4326", always_xy=True)
    lon, lat = to_ll.transform(xx.ravel(), yy.ravel())
    return np.asarray(lon).reshape(xx.shape), np.asarray(lat).reshape(yy.shape)


def remove_restore_mask(
    grid: GridLike,
    crs_epsg: int,
    areas: AreaDB | None = None,
    events: EventDB | None = None,
    reference_date: datetime | None = None,
    secondary_date: datetime | None = None,
) -> np.ndarray:
    """Boolean mask of pixels inside the areas to exclude (True = exclude).

    Parameters
    ----------
    grid : GridLike
        Target grid (paths, dataset or ``(x, y)`` coordinate arrays).
    crs_epsg : int
        EPSG of the grid coordinates (the GeoJSON is EPSG:4326).
    areas : AreaDB, optional
        Persistent areas; always applied.
    events : EventDB, optional
        Events; applied only to those that affect the pair, which needs
        both dates.
    reference_date, secondary_date : datetime, optional
        The pair's dates; required when `events` is given.

    """
    from shapely import contains_xy, union_all

    shapes: list[Any] = []
    if areas is not None:
        shapes.extend(areas.shapes())
    if events is not None and events.items:
        if reference_date is None or secondary_date is None:
            msg = "events need reference_date and secondary_date to decide which apply"
            raise ValueError(msg)
        active = events.active(reference_date, secondary_date)
        shapes.extend(e.shape() for e in active)
        logger.info("%d of %d events affect the pair", len(active), len(events.items))
    gx, gy = grid_coordinates(grid)
    shape_out = (len(np.asarray(gy)), len(np.asarray(gx)))
    if not shapes:
        return np.zeros(shape_out, dtype=bool)
    geom = union_all(shapes)
    lon, lat = _pixel_centres_lonlat(grid, crs_epsg)
    inside = contains_xy(geom, lon.ravel(), lat.ravel())
    mask = np.asarray(inside, dtype=bool).reshape(shape_out)
    logger.info(
        "remove-restore: %d pixels (%.2f %%) excluded from the fit",
        int(mask.sum()),
        100 * mask.mean(),
    )
    return mask


def sigma_inflation_inside(
    mask: np.ndarray,
    scale_px: float,
    max_factor: float = 3.0,
) -> np.ndarray:
    """Factor by which sigma grows inside excluded areas (PRD R-E1).

    1 outside and on the boundary, rising with the distance into the area as
    ``1 + (max_factor - 1) * (1 - exp(-d / scale_px))``: deep inside a basin
    the surface is an interpolation across it and its uncertainty saturates
    at `max_factor` times the fitted sigma.
    """
    from scipy.ndimage import distance_transform_edt

    m = np.asarray(mask, dtype=bool)
    if scale_px <= 0 or max_factor < 1:
        msg = "scale_px must be > 0 and max_factor >= 1"
        raise ValueError(msg)
    if not m.any():
        return np.ones(m.shape, dtype=np.float64)
    d = distance_transform_edt(m)
    return 1.0 + (max_factor - 1.0) * (1.0 - np.exp(-d / scale_px))


@dataclass
class RemoveRestore:
    """Everything the fit needs from the curated databases, for one pair."""

    mask: np.ndarray
    area_version: str | None = None
    event_version: str | None = None
    active_event_ids: list[str] = field(default_factory=list)

    @classmethod
    def for_pair(
        cls,
        grid: GridLike,
        crs_epsg: int,
        reference_date: datetime,
        secondary_date: datetime,
        areas: AreaDB | None = None,
        events: EventDB | None = None,
    ) -> RemoveRestore:
        mask = remove_restore_mask(
            grid, crs_epsg, areas, events, reference_date, secondary_date
        )
        active = events.active(reference_date, secondary_date) if events else []
        return cls(
            mask=mask,
            area_version=areas.version if areas else None,
            event_version=events.version if events else None,
            active_event_ids=[e.id for e in active],
        )
