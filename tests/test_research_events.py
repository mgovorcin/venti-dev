# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Event footprints and the bundled area / event databases (plan T41)."""

from __future__ import annotations

from datetime import datetime
from importlib import resources

import pytest
from shapely.geometry import Point, shape

from venti.calibration.remove_restore import load_area_db, load_event_db
from venti.research.events import event_feature, rupture_length_km

DATA = resources.files("venti") / "data"


def test_rupture_length_scaling():
    assert rupture_length_km(7.1) == pytest.approx(50.6, abs=0.2)  # Ridgecrest ~50 km
    assert rupture_length_km(6.0) < rupture_length_km(6.4) < rupture_length_km(7.0)


def test_footprint_contains_the_epicentre_and_follows_the_strike():
    f = event_feature(
        "e",
        "test",
        datetime(2020, 1, 1),
        -117.0,
        35.0,
        7.0,
        0.0,
        source="test",
        min_buffer_km=5.0,
        buffer_factor=0.1,
    )
    poly = shape(f["geometry"])
    assert poly.is_valid
    assert poly.contains(Point(-117.0, 35.0))
    minx, miny, maxx, maxy = poly.bounds
    assert (maxy - miny) > 3 * (maxx - minx) * 0.82  # elongated N-S (strike 0)
    assert f["properties"]["buffer_km"] == pytest.approx(5.0)


@pytest.mark.parametrize(
    ("name", "loader"),
    [
        ("defo_area_db_v1.geojson", load_area_db),
        ("event_db_v1.geojson", load_event_db),
    ],
)
def test_bundled_databases_are_valid(name, loader):
    db = loader(DATA / name)
    assert db.version
    ids = [i.id for i in db.items]
    assert ids
    assert len(ids) == len(set(ids))
    for item in db.items:
        assert item.shape().is_valid
        assert not item.shape().is_empty
        assert item.source


def test_ridgecrest_events_affect_only_spanning_pairs():
    db = load_event_db(DATA / "event_db_v1.geojson")
    span = db.active(datetime(2019, 6, 30), datetime(2019, 7, 12))
    after = db.active(datetime(2019, 7, 12), datetime(2019, 7, 24))
    assert {e.id for e in span} == {"ridgecrest_2019_mw71", "ridgecrest_2019_mw64"}
    assert after == []
    mw71 = next(e for e in db.items if e.id == "ridgecrest_2019_mw71")
    assert mw71.shape().contains(Point(-117.5993333, 35.7695))


def test_area_db_covers_the_benchmark_areas_and_reaches_stable_ground():
    """T41.2/T41.4: Houston, San Joaquin, Kern and Kilauea are in the DB, and
    every drafted polygon's 10 km ring is mostly stable ground."""
    db = load_area_db(DATA / "defo_area_db_v1.geojson")
    ids = {a.id for a in db.items}
    assert {
        "houston_galveston_subsidence",
        "san_joaquin_subsidence",
        "kern_subsidence_north",
        "kern_subsidence_south",
        "kilauea",
    } <= ids
    drafted = [a for a in db.items if hasattr(a, "stable_ring_fraction")]
    assert len(drafted) == 4
    assert all(a.stable_ring_fraction >= 0.5 for a in drafted)
    kilauea = next(a for a in db.items if a.id == "kilauea")
    assert kilauea.shape().contains(Point(-155.28, 19.41))  # the summit caldera
