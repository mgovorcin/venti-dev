# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Build a per-frame stable-pixel fit mask (plan T62, Houston/Matagorda).

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
"""

from __future__ import annotations

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np

__all__ = ["mean_coherence", "stable_pixel_mask", "write_fit_mask"]

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
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    products = sorted(a.products)
    mask = stable_pixel_mask(mean_coherence(products, a.layer), a.threshold)
    write_fit_mask(a.out, mask, products[0])
    logger.info(
        "%s: %d products, %.1f%% of the frame kept",
        a.out, len(products), 100 * mask.mean(),
    )  # fmt: skip


if __name__ == "__main__":  # pragma: no cover
    main()
