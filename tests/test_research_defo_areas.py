# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Draft deforming-area polygons (plan T41.1)."""

from __future__ import annotations

import json

import numpy as np
from pyproj import Transformer
from rasterio.transform import from_origin
from shapely.geometry import Point, shape

from venti.calibration.remove_restore import AreaDB
from venti.research.defo_areas import DraftOptions, draft_areas, main, regional_plane

X0, Y0, PIX = 500_000.0, 4_000_000.0, 200.0


def _scene(ny=400, nx=500):
    """A 30 mm/yr subsidence bowl (sigma 4 km) on a tilted plane plus noise."""
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:ny, 0:nx]
    plane = 3.0 + 0.01 * xx - 0.005 * yy
    cy, cx = 150, 300
    bowl = -30.0 * np.exp(
        -(((xx - cx) ** 2 + (yy - cy) ** 2) * PIX**2) / (2 * 4000.0**2)
    )
    v = plane + bowl + rng.normal(0, 0.5, (ny, nx))
    v[:20, :20] = np.nan
    x_c, y_c = X0 + (cx + 0.5) * PIX, Y0 - (cy + 0.5) * PIX
    return v, from_origin(X0, Y0, PIX, PIX), (x_c, y_c), plane


def test_plane_is_removed_robustly():
    v, _, _, plane = _scene()
    fit = regional_plane(v)
    far = np.isfinite(v)
    far[100:200, 250:350] = False  # the bowl
    assert np.nanmedian(np.abs(fit - plane)[far]) < 0.2


def test_one_buffered_polygon_around_the_bowl_reaching_stable_ground():
    v, tr, (xc, yc), _ = _scene()
    fc = draft_areas(
        v,
        tr,
        "EPSG:32611",
        DraftOptions(
            threshold_mm_yr=5,
            buffer_m=5000,
            min_area_km2=5,
            downsample=1,
            version="test-1",
            name_prefix="bowl",
        ),
    )
    assert fc["version"] == "test-1"
    assert len(fc["features"]) == 1
    props = fc["features"][0]["properties"]
    assert props["id"] == "bowl_01"
    assert props["reaches_stable_ground"]
    assert props["max_abs_residual_mm_yr"] > 20
    lon, lat = Transformer.from_crs(
        "EPSG:32611", "EPSG:4326", always_xy=True
    ).transform(xc, yc)
    poly = shape(fc["features"][0]["geometry"])
    assert poly.contains(Point(lon, lat))
    db = AreaDB.from_geojson(fc)  # the remove-restore schema reads it
    assert db.version == "test-1"
    assert db.items[0].id == "bowl_01"


def test_flat_field_gives_no_area_and_cli_writes_geojson(tmp_path):
    import rasterio

    rng = np.random.default_rng(1)
    v = 2.0 + rng.normal(0, 0.5, (200, 200))
    path = tmp_path / "v.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=200,
        width=200,
        count=1,
        dtype="float32",
        crs="EPSG:32611",
        transform=from_origin(X0, Y0, PIX, PIX),
    ) as ds:
        ds.write(v.astype(np.float32), 1)
    out = tmp_path / "draft.geojson"
    main([str(path), "--out", str(out), "--downsample", "1"])
    fc = json.loads(out.read_text())
    assert fc["features"] == []
    assert fc["version"] == "draft"
