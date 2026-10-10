# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Draft deforming-area polygons from a velocity raster (plan T41.1).

The remove-restore step (`venti.calibration.remove_restore`) needs curated,
versioned GeoJSON polygons of persistent deformation (basins, volcanoes).
This tool drafts them for a human to review, following the trade-study
recipe (rr4):

1. remove a regional trend from the LOS velocity (a robust plane); use a
   **calibrated** velocity: raw DISP keeps long-wavelength errors a plane
   cannot remove;
2. keep pixels whose residual exceeds `threshold_mm_yr` in magnitude,
   clean them with a morphological opening, keep connected areas of at least
   `min_area_km2`;
3. buffer each area by `buffer_m` (the polygon must reach stable ground, so
   the fit sees the ground around it), merge overlaps, simplify;
4. check that the buffer ring around each polygon is mostly stable
   (|residual| < threshold / 2 on at least `min_stable_fraction` of it) and
   record the result; a polygon that fails needs a larger buffer or a
   manual edit.

Output: a GeoJSON FeatureCollection in EPSG:4326 with a top-level
``version`` and per-feature ``id``, ``name``, ``source`` (the schema of
`venti.calibration.remove_restore.AreaDB`), plus review properties.

    venti-defo-areas VELOCITY.tif --out draft.geojson --threshold 5 --buffer 10000
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

__all__ = ["DraftOptions", "draft_areas", "main", "regional_plane"]


@dataclass
class DraftOptions:
    """Thresholds of the draft (defaults: the trade-study rr4 values)."""

    threshold_mm_yr: float = 5.0
    min_area_km2: float = 25.0
    buffer_m: float = 10_000.0
    simplify_m: float = 500.0
    opening_px: int = 3
    min_stable_fraction: float = 0.5
    downsample: int = 4
    sign: str = "both"
    version: str = "draft"
    name_prefix: str = "area"
    source: str = "venti-defo-areas draft from an InSAR velocity raster"


def regional_plane(v: np.ndarray, rng_seed: int = 0) -> np.ndarray:
    """Robust plane fit of `v` (two passes, 3-sigma clip), evaluated everywhere."""
    ny, nx = v.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    ok = np.isfinite(v)
    idx = np.flatnonzero(ok)
    rng = np.random.default_rng(rng_seed)
    if idx.size > 200_000:
        idx = rng.choice(idx, 200_000, replace=False)
    a = np.c_[np.ones(idx.size), xx.ravel()[idx], yy.ravel()[idx]]
    b = v.ravel()[idx]
    keep = np.ones(idx.size, dtype=bool)
    for _ in range(2):
        coef, *_ = np.linalg.lstsq(a[keep], b[keep], rcond=None)
        r = b - a @ coef
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        keep = np.abs(r) < 3 * s
    return coef[0] + coef[1] * xx + coef[2] * yy


def draft_areas(
    velocity: np.ndarray,
    transform: Any,
    crs: Any,
    options: DraftOptions | None = None,
) -> dict[str, Any]:
    """Draft area polygons from a LOS velocity grid (mm/yr, NaN = no data).

    Returns a GeoJSON FeatureCollection (EPSG:4326) for review.
    """
    from pyproj import Transformer
    from rasterio.features import shapes
    from scipy import ndimage
    from shapely.geometry import mapping, shape
    from shapely.ops import transform as shp_transform
    from shapely.ops import unary_union

    o = options or DraftOptions()
    f = max(1, int(o.downsample))
    v = np.asarray(velocity, dtype=np.float64)[::f, ::f]
    tr = transform * transform.scale(f, f) if f > 1 else transform
    pixel_m = abs(tr.a)
    resid = v - regional_plane(v)
    if o.sign == "negative":  # subsidence only
        anomaly = resid < -o.threshold_mm_yr
    elif o.sign == "positive":  # uplift only
        anomaly = resid > o.threshold_mm_yr
    else:
        anomaly = np.abs(resid) > o.threshold_mm_yr
    anomaly &= np.isfinite(resid)
    if o.opening_px > 0:
        anomaly = ndimage.binary_opening(anomaly, iterations=o.opening_px)
    lab, n = ndimage.label(anomaly, structure=np.ones((3, 3), bool))
    area_km2 = np.bincount(lab.ravel(), minlength=n + 1) * (pixel_m / 1000) ** 2
    big = [k for k in range(1, n + 1) if area_km2[k] >= o.min_area_km2]
    keep = np.isin(lab, big)

    polys = [
        shape(geom)
        for geom, val in shapes(keep.astype(np.uint8), mask=keep, transform=tr)
        if val == 1
    ]
    if not polys:
        logger.info(
            "no anomaly above %.1f mm/yr and %.0f km2",
            o.threshold_mm_yr,
            o.min_area_km2,
        )
        merged = []
    else:
        merged_geom = unary_union([p.buffer(o.buffer_m) for p in polys]).simplify(
            o.simplify_m
        )
        merged = list(getattr(merged_geom, "geoms", [merged_geom]))

    to_ll = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform
    stable = np.abs(resid) < o.threshold_mm_yr / 2
    finite = np.isfinite(resid)
    features = []
    for i, poly in enumerate(sorted(merged, key=lambda p: -p.area), 1):
        ring = poly.buffer(o.buffer_m).difference(poly)
        ring_mask = _rasterize(ring, v.shape, tr) & finite
        frac = float(stable[ring_mask].mean()) if ring_mask.any() else 0.0
        inside = _rasterize(poly, v.shape, tr) & finite
        features.append(
            {
                "type": "Feature",
                "properties": {
                    "id": f"{o.name_prefix}_{i:02d}",
                    "name": f"{o.name_prefix} {i}",
                    "source": o.source,
                    "area_km2": round(poly.area / 1e6, 1),
                    "max_abs_residual_mm_yr": round(
                        (
                            float(np.nanmax(np.abs(resid[inside])))
                            if inside.any()
                            else 0.0
                        ),
                        1,
                    ),
                    "stable_ring_fraction": round(frac, 2),
                    "reaches_stable_ground": frac >= o.min_stable_fraction,
                    "review": "draft: check extent, GNSS stations, buffer",
                },
                "geometry": mapping(shp_transform(to_ll, poly)),
            }
        )
    return {
        "type": "FeatureCollection",
        "version": o.version,
        "properties": {
            "threshold_mm_yr": o.threshold_mm_yr,
            "buffer_m": o.buffer_m,
            "min_area_km2": o.min_area_km2,
        },
        "features": features,
    }


def _rasterize(geom: Any, shape: tuple[int, int], transform: Any) -> np.ndarray:
    from rasterio.features import rasterize

    if geom.is_empty:
        return np.zeros(shape, dtype=bool)
    return rasterize([(geom, 1)], out_shape=shape, transform=transform, fill=0).astype(
        bool
    )


def main(argv: list[str] | None = None) -> None:
    """Command line: velocity GeoTIFF in, draft GeoJSON out."""
    import argparse

    import rasterio

    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("velocity", type=Path, help="LOS velocity GeoTIFF (mm/yr)")
    p.add_argument("--out", type=Path, required=True)
    d = DraftOptions()
    p.add_argument("--threshold", type=float, default=d.threshold_mm_yr)
    p.add_argument("--min-area", type=float, default=d.min_area_km2, help="km2")
    p.add_argument("--buffer", type=float, default=d.buffer_m, help="m")
    p.add_argument("--simplify", type=float, default=d.simplify_m, help="m")
    p.add_argument("--downsample", type=int, default=d.downsample)
    p.add_argument(
        "--sign",
        choices=("both", "negative", "positive"),
        default=d.sign,
        help="anomalies of either sign, subsidence only, or uplift only",
    )
    p.add_argument("--name", default=d.name_prefix, help="feature id/name prefix")
    p.add_argument("--version", default=d.version)
    a = p.parse_args(argv)
    with rasterio.open(a.velocity) as ds:
        v = ds.read(1, masked=True).filled(np.nan).astype(np.float64)
        transform, crs = ds.transform, ds.crs
    fc = draft_areas(
        v,
        transform,
        crs,
        DraftOptions(
            threshold_mm_yr=a.threshold,
            min_area_km2=a.min_area,
            buffer_m=a.buffer,
            simplify_m=a.simplify,
            downsample=a.downsample,
            sign=a.sign,
            name_prefix=a.name,
            version=a.version,
        ),
    )
    a.out.write_text(json.dumps(fc, indent=1))
    for feat in fc["features"]:
        pr = feat["properties"]
        print(
            f"{pr['id']}: {pr['area_km2']} km2, max |residual| "
            f"{pr['max_abs_residual_mm_yr']} mm/yr, stable ring "
            f"{pr['stable_ring_fraction']:.0%}"
            + (
                ""
                if pr["reaches_stable_ground"]
                else "  <- does NOT reach stable ground"
            )
        )


if __name__ == "__main__":
    main()
