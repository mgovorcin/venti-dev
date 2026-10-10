# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Event footprints for the event database (plan T41.3).

An earthquake's footprint is where its coseismic displacement is large
enough to bias a calibration surface fitted across it. The draft rule: the
surface rupture (length from the magnitude, Wells & Coppersmith 1994,
strike-slip SRL: ``log10 L_km = -3.55 + 0.74 M``,
https://doi.org/10.1785/BSSA0840040974) centred on the epicentre along the
strike, buffered by ``buffer_factor * L`` (at least `min_buffer_km`). A human
checks it against the coseismic step in a quick-look velocity (T41.4).

An eruption's footprint (`eruption_feature`) is the line through its vents
buffered by a fixed distance: the dike opening and flank motion follow the
rift, not a magnitude scaling. It has a ``t1``, so every pair overlapping
the eruption (not only those spanning one instant) excludes it.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

__all__ = ["eruption_feature", "event_feature", "rupture_length_km"]


def rupture_length_km(magnitude: float) -> float:
    """Surface rupture length of a strike-slip event (Wells & Coppersmith 1994)."""
    return 10 ** (-3.55 + 0.74 * magnitude)


def event_feature(
    event_id: str,
    name: str,
    t0: datetime,
    lon: float,
    lat: float,
    magnitude: float,
    strike_deg: float,
    *,
    source: str,
    buffer_factor: float = 0.75,
    min_buffer_km: float = 15.0,
    t1: datetime | None = None,
) -> dict[str, Any]:
    """Return the GeoJSON Feature of one event (`remove_restore.Event` schema)."""
    from pyproj import Geod
    from shapely.geometry import LineString, mapping
    from shapely.ops import transform

    length = rupture_length_km(magnitude)
    buffer_km = max(min_buffer_km, buffer_factor * length)
    geod = Geod(ellps="WGS84")
    ends = [
        geod.fwd(lon, lat, strike_deg + az, 500.0 * length)[:2] for az in (0.0, 180.0)
    ]
    # buffer in a local azimuthal-equidistant frame, metres
    from pyproj import Transformer

    aeqd = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +units=m +datum=WGS84"
    fwd = Transformer.from_crs("EPSG:4326", aeqd, always_xy=True).transform
    inv = Transformer.from_crs(aeqd, "EPSG:4326", always_xy=True).transform
    footprint = transform(fwd, LineString(ends)).buffer(
        buffer_km * 1000.0, quad_segs=16
    )
    return {
        "type": "Feature",
        "properties": {
            "id": event_id,
            "name": name,
            "source": source,
            "t0": t0.isoformat(),
            "t1": t1.isoformat() if t1 else None,
            "magnitude": magnitude,
            "epicentre": [lon, lat],
            "strike_deg": strike_deg,
            "rupture_length_km": round(length, 1),
            "buffer_km": round(buffer_km, 1),
            "footprint_rule": (
                "rupture (Wells & Coppersmith 1994 SRL) along strike through the"
                f" epicentre, buffered by max({min_buffer_km:g} km,"
                f" {buffer_factor:g} L)"
            ),
        },
        "geometry": mapping(transform(inv, footprint)),
    }


def eruption_feature(
    event_id: str,
    name: str,
    t0: datetime,
    t1: datetime,
    vents: list[tuple[float, float]],
    *,
    buffer_km: float,
    source: str,
) -> dict[str, Any]:
    """Return the GeoJSON Feature of an eruption: `vents` (lon, lat) buffered.

    The buffer is applied in a local azimuthal-equidistant frame centred on
    the vents' mean position, in kilometres.
    """
    from pyproj import Transformer
    from shapely.geometry import LineString, Point, mapping
    from shapely.ops import transform

    if not vents:
        msg = "an eruption needs at least one vent"
        raise ValueError(msg)
    if t1 < t0:
        msg = f"t1 {t1} is before t0 {t0}"
        raise ValueError(msg)
    lon0 = sum(v[0] for v in vents) / len(vents)
    lat0 = sum(v[1] for v in vents) / len(vents)
    aeqd = f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m +datum=WGS84"
    fwd = Transformer.from_crs("EPSG:4326", aeqd, always_xy=True).transform
    inv = Transformer.from_crs(aeqd, "EPSG:4326", always_xy=True).transform
    trace = LineString(vents) if len(vents) > 1 else Point(vents[0])
    footprint = transform(fwd, trace).buffer(buffer_km * 1000.0, quad_segs=16)
    return {
        "type": "Feature",
        "properties": {
            "id": event_id,
            "name": name,
            "source": source,
            "t0": t0.isoformat(),
            "t1": t1.isoformat(),
            "vents": [list(v) for v in vents],
            "buffer_km": buffer_km,
            "footprint_rule": f"line through the vents, buffered by {buffer_km:g} km",
        },
        "geometry": mapping(transform(inv, footprint)),
    }
