# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Water-mask regions for the unwrap-error correction (PRD R-U1; plan T36.1).

An unwrapping error is constant over a piece of land the unwrapper could not
connect to the mainland: an island or a strip cut off by a bay. The region
definition is therefore the *water mask*, not the displacement: erode the
land by three pixels so slivers and spits too thin to carry a measurement
vanish, grow the surviving seeds back to the coastline with a watershed on
the land mask, drop pieces under `min_region_area` pixels and label the
connected land that is left (the "Brisbane recipe"; 33 regions on F08882
instead of thousands of recommended-mask islands or 32 connected components
with 42 % of the pixels unlabelled). Land joined by a neck stays one region:
the unwrapper connects across coherent land, however narrow.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy import ndimage
from skimage import measure, morphology, segmentation

logger = logging.getLogger(__name__)

__all__ = ["downsample_labels", "largest_region", "segment_regions"]


def segment_regions(
    water_mask: np.ndarray,
    valid_mask: np.ndarray | None = None,
    *,
    erosion_iterations: int = 3,
    min_region_area: int = 20,
) -> np.ndarray:
    """Label land regions separated by water (0 = water or dropped).

    Parameters
    ----------
    water_mask : np.ndarray
        True on land (the DISP ``water_mask`` layer convention).
    valid_mask : np.ndarray, optional
        Pixels that carry data (e.g. finite ``temporal_coherence``); land
        without data is treated as water so no-data strips do not join
        regions.
    erosion_iterations : int
        Erosions by a 3x3 disk before seeding; pieces thinner than about
        twice this many pixels leave no seed and are dropped.
    min_region_area : int
        Regions smaller than this many pixels are dropped.

    Returns
    -------
    np.ndarray
        int32 labels, consecutive from 1, 0 elsewhere.

    """
    land = np.asarray(water_mask, dtype=bool)
    if land.ndim != 2:
        msg = f"water_mask must be 2-D, got shape {land.shape}"
        raise ValueError(msg)
    if valid_mask is not None:
        valid = np.asarray(valid_mask, dtype=bool)
        if valid.shape != land.shape:
            msg = f"valid_mask {valid.shape} does not match water_mask {land.shape}"
            raise ValueError(msg)
        land = land & valid
    if erosion_iterations < 0 or min_region_area < 1:
        msg = "erosion_iterations must be >= 0 and min_region_area >= 1"
        raise ValueError(msg)
    eroded = land.copy()
    disk = morphology.disk(1)
    for _ in range(erosion_iterations):
        eroded = ndimage.binary_erosion(eroded, disk)
    seeds = measure.label(eroded)
    ws = segmentation.watershed(-land.astype(float), seeds, mask=land)
    ids, counts = np.unique(ws[ws > 0], return_counts=True)
    ws[np.isin(ws, ids[counts < min_region_area])] = 0
    labels = measure.label(ws > 0).astype(np.int32)
    logger.info(
        "segment_regions: %d regions from %d seeds (%d dropped as < %d px)",
        int(labels.max()),
        int(seeds.max()),
        int((counts < min_region_area).sum()),
        min_region_area,
    )
    return labels


def largest_region(labels: np.ndarray) -> int:
    """Return the label of the largest region (the mainland), 0 if none."""
    ids, counts = np.unique(labels[labels > 0], return_counts=True)
    return int(ids[np.argmax(counts)]) if ids.size else 0


def downsample_labels(labels: np.ndarray, factor: int) -> np.ndarray:
    """Pick the block-centre label so labels follow `downsample_array`'s grid.

    The output shape is ``(ny // factor, nx // factor)``, the trimmed shape
    `venti.spatial.resample.downsample_array` produces.
    """
    lab = np.asarray(labels)
    if factor < 1:
        msg = "factor must be >= 1"
        raise ValueError(msg)
    if factor == 1:
        return lab
    ny, nx = lab.shape[0] // factor, lab.shape[1] // factor
    half = factor // 2
    return lab[half::factor, half::factor][:ny, :nx]
