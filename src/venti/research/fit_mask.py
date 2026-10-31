# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
r"""Build a per-frame stable-pixel fit mask (plan T62, Houston/Matagorda).

The loclin surface absorbs the phase-bias signal of decorrelating fields
(crops, wetlands) around stable GNSS towns. A frame whose benchmark shows it
gets a fixed mask of pixels that are coherent over the whole archive; the
fit uses only those pixels (`SurfaceOptions.fit_mask`), and the surface is
still evaluated everywhere.

Recipe: the mean of every DISP product's ``temporal_coherence`` (NaN
skipped), thresholded. The stack mean matters: a per-pair threshold keeps
fields that are coherent for a while (20-30% of the Houston frame at 0.7 in
early pairs, against 9% for the stack mean) and scored worse in the T46 edge
pass.

    venti-fit-mask DISP_DIR/*.nc --threshold 0.7 --out F08882_stable_coh07.tif

A high threshold removes more of that signal but leaves little in rural parts
of the frame, where the surface is then extrapolated. `grow_fit_mask` keeps
both: seeds at a high threshold, then, at falling coherence, pixels close to
the mask whose stack velocity matches the seeds' local velocity and whose mean
phase similarity is high. The inputs are DISP only; no GNSS.

    venti-fit-mask DISP_DIR/*.nc --threshold 0.85 --grow-floor 0.6 \\
        --velocity VELOCITY.tif --water-mask WATER.tif --area-db AREAS.geojson \\
        --out F38238_stable_g085_06.tif
"""

from __future__ import annotations

import argparse
import logging
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np

__all__ = ["grow_fit_mask", "mean_coherence", "stable_pixel_mask", "write_fit_mask"]

logger = logging.getLogger(__name__)


def mean_coherence(
    paths: Sequence[Path | str], layer: str = "temporal_coherence"
) -> np.ndarray:
    """Mean of `layer` over the DISP products in `paths` (NaN where never valid)."""
    import h5py

    total = count = None
    for p in paths:
        with h5py.File(p, "r") as f:
            c = f[layer][()].astype(np.float32)
        ok = np.isfinite(c)
        if total is None:
            total = np.zeros(c.shape, np.float64)
            count = np.zeros(c.shape, np.int32)
        elif c.shape != total.shape:
            msg = f"{Path(p).name}: {layer} shape {c.shape} != {total.shape}"
            raise ValueError(msg)
        total[ok] += c[ok]
        count += ok
    if total is None or count is None:
        msg = "no products given"
        raise ValueError(msg)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(count > 0, total / count, np.nan).astype(np.float32)


def stable_pixel_mask(coherence: np.ndarray, threshold: float) -> np.ndarray:
    """uint8 mask, 1 where the mean coherence is at least `threshold`."""
    return (np.nan_to_num(coherence, nan=0.0) >= threshold).astype(np.uint8)


def _local_seed_velocity(
    velocity: np.ndarray,
    seeds: np.ndarray,
    block_px: int,
    sigma_blocks: float,
) -> np.ndarray:
    """Smoothed median velocity of the seeds, on the grid of `velocity`.

    Block medians are spread by a normalized Gaussian, NaN where no seed is
    within reach. Pixels past the last whole block take its value.
    """
    from scipy.ndimage import gaussian_filter

    ny, nx = velocity.shape
    h, w = ny // block_px * block_px, nx // block_px * block_px
    blocks = np.where(seeds, velocity, np.nan)[:h, :w].reshape(
        h // block_px, block_px, w // block_px, block_px
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # blocks without a seed
        median = np.nanmedian(blocks, axis=(1, 3))
    ok = np.isfinite(median)
    weight = gaussian_filter(ok.astype(float), sigma_blocks)
    smooth = gaussian_filter(np.where(ok, median, 0.0), sigma_blocks)
    smooth = smooth / np.maximum(weight, 1e-3)
    smooth[weight < 1e-3] = np.nan
    rows = np.minimum(np.arange(ny), h - 1) // block_px
    cols = np.minimum(np.arange(nx), w - 1) // block_px
    return smooth[np.ix_(rows, cols)]


def grow_fit_mask(
    coherence: np.ndarray,
    velocity: np.ndarray,
    *,
    seed_threshold: float = 0.85,
    floor: float = 0.6,
    step: float = 0.05,
    similarity: np.ndarray | None = None,
    min_similarity: float = 0.8,
    max_velocity_diff: float = 2.0,
    reach_m: float = 2000.0,
    valid: np.ndarray | None = None,
    exclude: np.ndarray | None = None,
    posting_m: float = 30.0,
    velocity_block_m: float = 1000.0,
    velocity_sigma_m: float = 5000.0,
) -> np.ndarray:
    """Grow a fit mask from high-coherence seeds.

    Parameters
    ----------
    coherence : np.ndarray
        Stack-mean temporal coherence (`mean_coherence`).
    velocity : np.ndarray
        Stack LOS velocity, in the units of `max_velocity_diff`; NaN where
        invalid. Only local differences matter, so its reference does not.
    seed_threshold, floor, step : float
        Seeds have coherence >= `seed_threshold`; candidates are added at
        ``seed_threshold - step``, then lower, down to `floor`.
    similarity : np.ndarray, optional
        Stack-mean ``phase_similarity``; candidates need `min_similarity`.
    min_similarity : float
        Used with `similarity`.
    max_velocity_diff : float
        A candidate's velocity must be this close to the seeds' local velocity.
    reach_m : float
        A candidate must be this close to the mask of the previous level.
    valid : np.ndarray, optional
        True where a pixel may be used (e.g. land); pixels with NaN velocity
        never are.
    exclude : np.ndarray, optional
        True on known deforming areas. They stay in the mask if they qualify
        but do not set the local velocity.
    posting_m : float
        Pixel size.
    velocity_block_m, velocity_sigma_m : float
        Block size and Gaussian sigma of the local velocity.

    Returns
    -------
    np.ndarray
        uint8 mask, 1 = fit this pixel.

    """
    from scipy.ndimage import distance_transform_edt

    if coherence.shape != velocity.shape:
        msg = f"coherence {coherence.shape} != velocity {velocity.shape}"
        raise ValueError(msg)
    usable = np.isfinite(velocity)
    if valid is not None:
        usable &= valid
    coh = np.nan_to_num(coherence, nan=0.0)
    seeds = (coh >= seed_threshold) & usable
    ref_seeds = seeds if exclude is None else seeds & ~exclude
    block_px = int(velocity_block_m / posting_m)
    local = _local_seed_velocity(
        velocity, ref_seeds, block_px, velocity_sigma_m / velocity_block_m
    )
    with np.errstate(invalid="ignore"):
        candidate = usable & (np.abs(velocity - local) <= max_velocity_diff)
    if similarity is not None:
        candidate &= np.nan_to_num(similarity, nan=0.0) >= min_similarity

    mask = seeds.copy()
    factor = 4  # distances on a coarser grid: reach is far above a pixel
    for level in np.arange(seed_threshold - step, floor - 1e-6, -step):
        distance = (
            distance_transform_edt(~mask[::factor, ::factor]) * posting_m * factor
        )
        near = np.repeat(np.repeat(distance <= reach_m, factor, 0), factor, 1)
        near = near[: mask.shape[0], : mask.shape[1]]
        mask |= (coh >= level) & candidate & near
    return mask.astype(np.uint8)


def _read_on_grid(path: Path | str, shape: tuple[int, int]) -> np.ndarray:
    """Read band 1 of a raster that must already be on the mask grid."""
    import rasterio

    with rasterio.open(path) as src:
        data = src.read(1)
    if data.shape != shape:
        msg = f"{path}: shape {data.shape} != product grid {shape}"
        raise ValueError(msg)
    return data


def _area_raster(
    db_path: Path | str, product: Path | str, shape: tuple[int, int]
) -> np.ndarray:
    """Rasterize a Venti area database on the grid of the DISP `product`."""
    import h5py
    import shapely
    from affine import Affine
    from pyproj import CRS, Transformer
    from rasterio.features import rasterize

    from venti.calibration.remove_restore import load_area_db

    with h5py.File(product, "r") as f:
        x, y = f["x"][()], f["y"][()]
        wkt = f["spatial_ref"].attrs["crs_wkt"]
    wkt = wkt.decode() if isinstance(wkt, bytes) else str(wkt)
    to_grid = Transformer.from_crs(4326, CRS.from_wkt(wkt), always_xy=True).transform

    def project(lonlat: np.ndarray) -> np.ndarray:
        return np.column_stack(to_grid(lonlat[:, 0], lonlat[:, 1]))

    dx, dy = float(x[1] - x[0]), float(y[0] - y[1])
    return rasterize(
        [(shapely.transform(g, project), 1) for g in load_area_db(db_path).shapes()],
        out_shape=shape,
        transform=Affine(dx, 0.0, x[0] - dx / 2, 0.0, -dy, y[0] + dy / 2),
        fill=0,
        dtype="uint8",
    ).astype(bool)


def write_fit_mask(path: Path | str, mask: np.ndarray, product: Path | str) -> Path:
    """Write `mask` as a uint8 GeoTIFF on the grid of the DISP `product`."""
    import h5py
    import rasterio
    from affine import Affine

    with h5py.File(product, "r") as f:
        x, y = f["x"][()], f["y"][()]
        wkt = f["spatial_ref"].attrs["crs_wkt"]
    wkt = wkt.decode() if isinstance(wkt, bytes) else str(wkt)
    if mask.shape != (y.size, x.size):
        msg = f"mask shape {mask.shape} != product grid {(y.size, x.size)}"
        raise ValueError(msg)
    dx, dy = float(x[1] - x[0]), float(y[0] - y[1])
    transform = Affine(dx, 0.0, x[0] - dx / 2, 0.0, -dy, y[0] + dy / 2)
    path = Path(path)
    with rasterio.open(
        path, "w", driver="GTiff", height=mask.shape[0], width=mask.shape[1],
        count=1, dtype="uint8", crs=wkt, transform=transform,
        compress="deflate", predictor=2, tiled=True,
    ) as dst:  # fmt: skip
        dst.write(mask.astype(np.uint8), 1)
    return path


def main(argv: list[str] | None = None) -> None:
    """CLI: ``venti-fit-mask PRODUCTS... --threshold T --out MASK.tif``."""
    ap = argparse.ArgumentParser(prog="venti-fit-mask", description=__doc__)
    ap.add_argument("products", nargs="+", type=Path, help="DISP-S1 products")
    ap.add_argument("--threshold", type=float, default=0.7)
    ap.add_argument("--layer", default="temporal_coherence")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--grow-floor",
        type=float,
        help="grow from --threshold down to this coherence (needs --velocity)",
    )
    ap.add_argument(
        "--velocity", type=Path, help="stack LOS velocity (mm/yr), DISP grid"
    )
    ap.add_argument("--water-mask", type=Path, help="1 = land, DISP grid")
    ap.add_argument("--area-db", type=Path, help="Venti area GeoJSON")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    products = sorted(a.products)
    coherence = mean_coherence(products, a.layer)
    if a.grow_floor is None:
        mask = stable_pixel_mask(coherence, a.threshold)
    else:
        if a.velocity is None:
            ap.error("--grow-floor needs --velocity")
        grid = coherence.shape
        mask = grow_fit_mask(
            coherence,
            _read_on_grid(a.velocity, grid),
            seed_threshold=a.threshold,
            floor=a.grow_floor,
            similarity=mean_coherence(products, "phase_similarity"),
            valid=(
                None if a.water_mask is None else _read_on_grid(a.water_mask, grid) == 1
            ),
            exclude=(
                None
                if a.area_db is None
                else _area_raster(a.area_db, products[0], grid)
            ),
        )
    write_fit_mask(a.out, mask, products[0])
    logger.info(
        "%s: %d products, %.1f%% of the frame kept",
        a.out, len(products), 100 * mask.mean(),
    )  # fmt: skip


if __name__ == "__main__":  # pragma: no cover
    main()
