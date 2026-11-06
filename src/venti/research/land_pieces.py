# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Build a per-frame land-piece label map for breaklines (``surface.regions='pieces'``).

A per-pair calibration surface blurs a step between land pieces separated by
water (an unwrapping offset or disconnected components) up to ~20 km into
both (Gulf/East campaign: a 92 mm/yr LOS step in F28220). Each large piece of
a fixed map gets a surface of its own instead. The map is built once per
frame from the whole stack, so every pair of the frame uses the same pieces.

Recipe (the per-pair trade study, cal-vlm-batch gulf_east validate): pixels
that are land in the ``water_mask`` layer and have a stack-mean temporal
coherence of at least `--min-coherence`; blocks of `--factor` x `--factor`
pixels (the fit grid) more than half valid; gaps narrower than about twice
`--close-iterations` blocks bridged; connected pieces labelled
(`venti.unwrap.regions.segment_regions`). The labels are written at full
resolution on the frame's DISP grid.

    venti-land-pieces DISP_DIR/*.nc --out F28220_land_pieces.tif
"""

from __future__ import annotations

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from scipy import ndimage

from ..unwrap.regions import segment_regions
from .fit_mask import mean_coherence

__all__ = ["land_piece_labels", "land_pieces_from_products", "write_land_pieces"]

logger = logging.getLogger(__name__)


def land_piece_labels(
    valid: np.ndarray,
    *,
    factor: int = 6,
    close_iterations: int = 8,
    min_region_area: int = 50,
) -> np.ndarray:
    """Full-resolution piece labels (0 = no piece) from a full-resolution mask."""
    valid = np.asarray(valid, dtype=bool)
    h, w = valid.shape[0] // factor, valid.shape[1] // factor
    blocks = (
        valid[: h * factor, : w * factor]
        .reshape(h, factor, w, factor)
        .mean(axis=(1, 3))
        > 0.5
    )
    bridged = ndimage.binary_closing(
        blocks, np.ones((3, 3), bool), iterations=close_iterations
    )
    labels_ds = segment_regions(
        bridged, erosion_iterations=1, min_region_area=min_region_area
    )
    labels = np.zeros(valid.shape, dtype=np.uint16)
    labels[: h * factor, : w * factor] = np.repeat(
        np.repeat(labels_ds, factor, axis=0), factor, axis=1
    )
    return labels


def land_pieces_from_products(
    paths: Sequence[Path | str],
    *,
    min_coherence: float = 0.5,
    factor: int = 6,
    close_iterations: int = 8,
) -> np.ndarray:
    """Label map from the stack-mean coherence and the land mask of the DISP `paths`."""
    import h5py

    coherence = mean_coherence(paths)
    valid = np.nan_to_num(coherence, nan=0.0) >= min_coherence
    with h5py.File(paths[0], "r") as f:
        if "water_mask" in f:
            valid &= np.asarray(f["water_mask"][()], dtype=bool)
    return land_piece_labels(valid, factor=factor, close_iterations=close_iterations)


def write_land_pieces(
    path: Path | str, labels: np.ndarray, product: Path | str
) -> Path:
    """Write `labels` as a uint16 GeoTIFF on the grid of the DISP `product`."""
    import h5py
    import rasterio
    from affine import Affine

    with h5py.File(product, "r") as f:
        x, y = f["x"][()], f["y"][()]
        wkt = f["spatial_ref"].attrs["crs_wkt"]
    wkt = wkt.decode() if isinstance(wkt, bytes) else str(wkt)
    if labels.shape != (y.size, x.size):
        msg = f"labels shape {labels.shape} != product grid {(y.size, x.size)}"
        raise ValueError(msg)
    dx, dy = float(x[1] - x[0]), float(y[0] - y[1])
    transform = Affine(dx, 0.0, x[0] - dx / 2, 0.0, -dy, y[0] + dy / 2)
    path = Path(path)
    with rasterio.open(
        path, "w", driver="GTiff", height=labels.shape[0], width=labels.shape[1],
        count=1, dtype="uint16", nodata=0, crs=wkt, transform=transform,
        compress="deflate", predictor=2, tiled=True,
    ) as dst:  # fmt: skip
        dst.write(labels.astype(np.uint16), 1)
    return path


def main(argv: list[str] | None = None) -> None:
    """CLI: ``venti-land-pieces PRODUCTS... --out PIECES.tif``."""
    ap = argparse.ArgumentParser(prog="venti-land-pieces", description=__doc__)
    ap.add_argument("products", nargs="+", type=Path, help="DISP-S1 products")
    ap.add_argument("--min-coherence", type=float, default=0.5)
    ap.add_argument(
        "--factor",
        type=int,
        default=6,
        help="fit-grid block size (calibration downsample_factor)",
    )
    ap.add_argument("--close-iterations", type=int, default=8)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    products = sorted(a.products)
    labels = land_pieces_from_products(
        products,
        min_coherence=a.min_coherence,
        factor=a.factor,
        close_iterations=a.close_iterations,
    )
    write_land_pieces(a.out, labels, products[0])
    ids, counts = np.unique(labels[labels > 0], return_counts=True)
    sizes = sorted((counts // a.factor**2).tolist(), reverse=True)[:10]
    logger.info(
        "%s: %d products, %d pieces, sizes in fit cells %s",
        a.out, len(products), ids.size, sizes,
    )  # fmt: skip


if __name__ == "__main__":  # pragma: no cover
    main()
