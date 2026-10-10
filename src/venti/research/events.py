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
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

__all__ = ["event_feature", "rupture_length_km"]


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
