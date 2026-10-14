# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Stable-pixel fit mask recipe (venti-fit-mask, plan T62)."""

from __future__ import annotations

import h5py
import numpy as np
import pytest
import rasterio
from pyproj import CRS

from venti.research.fit_mask import main, mean_coherence, stable_pixel_mask


def _product(path, coh, x0=250_000.0, y0=3_300_000.0, posting=30.0):
    ny, nx = coh.shape
    with h5py.File(path, "w") as f:
        f["temporal_coherence"] = coh.astype(np.float32)
        f["x"] = x0 + posting * (np.arange(nx) + 0.5)
        f["y"] = y0 - posting * (np.arange(ny) + 0.5)
        ref = f.create_dataset("spatial_ref", data=0)
        ref.attrs["crs_wkt"] = CRS.from_epsg(32615).to_wkt()
    return path


def test_stack_mean_keeps_only_pixels_coherent_over_the_archive(tmp_path):
    """A field coherent in one pair of three (a crop between seasons) is not
    stable; NaN pairs are skipped and a never-valid pixel is not kept."""
    stack = [np.full((6, 8), 0.9) for _ in range(3)]
    stack[0][:, :4] = stack[1][:, :4] = 0.2  # crop field: 0.2, 0.2, 0.9
    stack[2][1, 6] = np.nan  # one missing pair: mean of the other two
    for c in stack:
        c[5, 7] = np.nan  # never valid
    paths = [_product(tmp_path / f"p{i}.nc", c) for i, c in enumerate(stack)]
    mean = mean_coherence(paths)
    assert mean[0, 0] == pytest.approx((0.2 + 0.2 + 0.9) / 3)
    assert mean[1, 6] == pytest.approx(0.9)
    assert np.isnan(mean[5, 7])
    mask = stable_pixel_mask(mean, 0.7)
    assert mask.dtype == np.uint8
    assert mask[:, :4].sum() == 0
    assert mask[:, 4:].sum() == 6 * 4 - 1
    assert mask[5, 7] == 0

    out = tmp_path / "mask.tif"
    main([*map(str, paths), "--threshold", "0.7", "--out", str(out)])
    with rasterio.open(out) as src:
        np.testing.assert_array_equal(src.read(1), mask)
        assert src.crs.to_epsg() == 32615
        assert src.transform.c == 250_000.0
        assert src.transform.f == 3_300_000.0
        assert src.transform.a == 30.0


def test_mean_coherence_rejects_mixed_grids(tmp_path):
    a = _product(tmp_path / "a.nc", np.ones((4, 4)))
    b = _product(tmp_path / "b.nc", np.ones((4, 5)))
    with pytest.raises(ValueError, match="shape"):
        mean_coherence([a, b])
