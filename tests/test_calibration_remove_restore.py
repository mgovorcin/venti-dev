# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Remove-restore areas and events (PRD A5, R-S5, §2.8; plan T32)."""

from __future__ import annotations

import json
from datetime import datetime

import numpy as np
import pytest
from pyproj import Transformer

from venti.calibration.gaps import base_weights, fill_gaps
from venti.calibration.loclin import loclin_surface
from venti.calibration.remove_restore import (
    AreaDB,
    EventDB,
    RemoveRestore,
    load_area_db,
    load_event_db,
    remove_restore_mask,
    sigma_inflation_inside,
)

UTM = 32615
X0, Y0, SIZE, N = 250_000.0, 3_300_000.0, 60_000.0, 120  # 500 m posting
PIXEL_M = SIZE / N


@pytest.fixture
def grid():
    x = X0 + (np.arange(N) + 0.5) * PIXEL_M
    y = Y0 + SIZE - (np.arange(N) + 0.5) * PIXEL_M
    return x, y


def box_lonlat(x0, y0, x1, y1):
    to_ll = Transformer.from_crs(f"EPSG:{UTM}", "EPSG:4326", always_xy=True)
    xs = [x0, x1, x1, x0, x0]
    ys = [y0, y0, y1, y1, y0]
    lon, lat = to_ll.transform(xs, ys)
    return {
        "type": "Polygon",
        "coordinates": [list(map(list, zip(lon, lat, strict=True)))],
    }


@pytest.fixture
def basin_geojson():
    return {
        "type": "FeatureCollection",
        "version": "defo-v1",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "id": "central_valley_test",
                    "name": "bowl",
                    "source": "test",
                },
                "geometry": box_lonlat(
                    X0 + 20_000, Y0 + 20_000, X0 + 40_000, Y0 + 40_000
                ),
            }
        ],
    }


@pytest.fixture
def events_geojson():
    return {
        "type": "FeatureCollection",
        "version": "events-v1",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "id": "ridgecrest_test",
                    "t0": "2019-07-06",
                    "magnitude": 7.1,
                    "source": "ComCat",
                },
                "geometry": box_lonlat(
                    X0 + 5_000, Y0 + 5_000, X0 + 15_000, Y0 + 15_000
                ),
            },
            {
                "type": "Feature",
                "properties": {
                    "id": "postseismic_test",
                    "t0": "2020-01-01",
                    "t1": "2020-12-31",
                },
                "geometry": box_lonlat(
                    X0 + 45_000, Y0 + 45_000, X0 + 55_000, Y0 + 55_000
                ),
            },
        ],
    }


def test_load_databases(tmp_path, basin_geojson, events_geojson):
    a = tmp_path / "defo.geojson"
    a.write_text(json.dumps(basin_geojson))
    e = tmp_path / "events.geojson"
    e.write_text(json.dumps(events_geojson))
    areas = load_area_db(a)
    events = load_event_db(e)
    assert areas.version == "defo-v1"
    assert [x.id for x in areas.items] == ["central_valley_test"]
    assert areas.shapes()[0].is_valid
    assert events.version == "events-v1"
    assert events.items[0].t0 == datetime(2019, 7, 6)
    assert events.items[1].t1 == datetime(2020, 12, 31)
    assert events.items[0].magnitude == 7.1


def test_version_is_required_and_type_checked():
    with pytest.raises(ValueError, match="version"):
        AreaDB.from_geojson({"type": "FeatureCollection", "features": []})
    with pytest.raises(ValueError, match="FeatureCollection"):
        EventDB.from_geojson({"type": "Feature", "version": "x"})


def test_event_spanning_rules(events_geojson):
    events = EventDB.from_geojson(events_geojson)
    coseismic, postseismic = events.items
    assert coseismic.affects(datetime(2019, 6, 1), datetime(2019, 8, 1))
    assert not coseismic.affects(datetime(2019, 8, 1), datetime(2019, 10, 1))  # after
    assert not coseismic.affects(datetime(2019, 1, 1), datetime(2019, 7, 5))  # before
    # a pair inside the postseismic window is affected even without spanning t0
    assert postseismic.affects(datetime(2020, 5, 1), datetime(2020, 6, 1))
    assert not postseismic.affects(datetime(2021, 2, 1), datetime(2021, 3, 1))
    assert [
        e.id for e in events.active(datetime(2019, 6, 1), datetime(2019, 8, 1))
    ] == ["ridgecrest_test"]


def test_mask_matches_the_polygons(grid, basin_geojson, events_geojson):
    areas = AreaDB.from_geojson(basin_geojson)
    events = EventDB.from_geojson(events_geojson)
    gx, gy = grid
    xx, yy = np.meshgrid(gx, gy)
    in_basin = (
        (xx > X0 + 20_000)
        & (xx < X0 + 40_000)
        & (yy > Y0 + 20_000)
        & (yy < Y0 + 40_000)
    )
    in_quake = (
        (xx > X0 + 5_000) & (xx < X0 + 15_000) & (yy > Y0 + 5_000) & (yy < Y0 + 15_000)
    )

    only_areas = remove_restore_mask(grid, UTM, areas)
    # the lon/lat polygon is the UTM box transformed, so agreement is to the pixel
    assert np.mean(only_areas == in_basin) > 0.995

    spanning = remove_restore_mask(
        grid, UTM, areas, events, datetime(2019, 6, 1), datetime(2019, 8, 1)
    )
    assert np.mean(spanning == (in_basin | in_quake)) > 0.995
    not_spanning = remove_restore_mask(
        grid, UTM, areas, events, datetime(2021, 1, 1), datetime(2021, 3, 1)
    )
    np.testing.assert_array_equal(
        not_spanning, only_areas
    )  # identical to "no event DB"
    assert not remove_restore_mask(grid, UTM).any()
    with pytest.raises(ValueError, match="reference_date"):
        remove_restore_mask(grid, UTM, areas, events)


def test_for_pair_records_versions_and_active_events(
    grid, basin_geojson, events_geojson
):
    rr = RemoveRestore.for_pair(
        grid,
        UTM,
        datetime(2019, 6, 1),
        datetime(2019, 8, 1),
        AreaDB.from_geojson(basin_geojson),
        EventDB.from_geojson(events_geojson),
    )
    assert rr.area_version == "defo-v1"
    assert rr.event_version == "events-v1"
    assert rr.active_event_ids == ["ridgecrest_test"]
    assert rr.mask.any()


def test_surface_is_interpolated_through_a_bowl(grid, basin_geojson):
    """A 20 cm bowl inside the area leaves the fitted surface flat there (R-S5)."""
    gx, gy = grid
    xx, yy = np.meshgrid(gx, gy)
    plane = 2.0 + 0.00002 * (xx - X0) - 0.00001 * (yy - Y0)  # mm, gentle regional field
    # a 3 km bowl: < 1 mm at the area boundary 10 km away, so the curated area
    # reaches stable ground as the PRD requires
    bowl = -200.0 * np.exp(
        -(((xx - X0 - 30_000) ** 2 + (yy - Y0 - 30_000) ** 2) / (2 * 3_000.0**2))
    )
    field = plane + bowl
    mask = remove_restore_mask(grid, UTM, AreaDB.from_geojson(basin_geojson))
    assert bowl[mask].min() < -150  # the bowl is where the area is

    valid = ~mask
    filled, filled_mask = fill_gaps(np.where(valid, field, np.nan), valid)
    w = base_weights(valid, filled_mask)
    surface, _ = loclin_surface(filled, w, PIXEL_M, 50_000.0)
    # inside the area the surface follows the regional plane, not the bowl
    assert np.abs(surface - plane)[mask].max() < 1.0
    # and the fit without exclusion would have swallowed part of it
    naive, _ = loclin_surface(field, np.ones_like(field), PIXEL_M, 50_000.0)
    assert np.abs(naive - plane)[mask].max() > 10.0


def test_sigma_inflation_inside():
    mask = np.zeros((40, 40), bool)
    mask[10:30, 10:30] = True
    f = sigma_inflation_inside(mask, scale_px=5.0, max_factor=3.0)
    assert f.shape == mask.shape
    np.testing.assert_array_equal(f[~mask], 1.0)
    assert f[20, 20] > f[10, 10] >= 1.0  # grows inward
    assert f.max() <= 3.0
    assert f[20, 20] == pytest.approx(1 + 2 * (1 - np.exp(-10 / 5)), rel=0.05)
    np.testing.assert_array_equal(
        sigma_inflation_inside(np.zeros((3, 3), bool), 5.0), 1.0
    )
    with pytest.raises(ValueError, match="scale_px"):
        sigma_inflation_inside(mask, scale_px=0)
