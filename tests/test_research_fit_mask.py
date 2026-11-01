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

from venti.research.fit_mask import (
    grow_fit_mask,
    main,
    mean_coherence,
    stable_pixel_mask,
)


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


def _grow_scene():
    """A 2 km block of seeds (coherence 0.9) in a 12 km scene at velocity 0."""
    n = 400
    coh = np.full((n, n), 0.4)
    vel = np.zeros((n, n))
    coh[180:220, 180:220] = 0.9
    return coh, vel


GROW = {"posting_m": 30.0, "reach_m": 300.0, "velocity_sigma_m": 2000.0}


def test_grow_adds_similar_pixels_near_the_seeds_level_by_level():
    coh, vel = _grow_scene()
    coh[220:226, 180:220] = 0.82  # touches the seeds: first level
    coh[226:232, 180:220] = 0.72  # reached through the first level only
    coh[240:246, 180:220] = 0.82  # coherent, but too far from level one's mask
    coh[300:306, 180:220] = 0.82  # far from everything
    mask = grow_fit_mask(coh, vel, seed_threshold=0.85, floor=0.7, **GROW)
    assert mask.dtype == np.uint8
    assert mask[180:220, 180:220].all()
    assert mask[220:226, 180:220].all()
    assert mask[226:232, 180:220].all()
    assert not mask[240:246, 180:220].any()
    assert not mask[300:306, 180:220].any()
    assert mask[coh == 0.4].sum() == 0


def test_grow_rejects_pixels_that_move_unlike_the_seeds():
    """Coherent subsiding ground stays out, whatever its coherence."""
    coh, vel = _grow_scene()
    coh[220:226, 180:200] = 0.8
    coh[220:226, 200:220] = 0.8
    vel[220:226, 200:220] = -10.0
    mask = grow_fit_mask(coh, vel, seed_threshold=0.85, floor=0.7, **GROW)
    assert mask[220:226, 180:200].all()
    assert not mask[220:226, 200:220].any()


def test_grow_honours_similarity_valid_and_nan():
    coh, vel = _grow_scene()
    coh[220:226, 180:220] = 0.8
    sim = np.ones_like(coh)
    sim[220:226, 180:190] = 0.5
    valid = np.ones_like(coh, dtype=bool)
    valid[220:226, 190:200] = False
    vel[220:226, 200:210] = np.nan
    mask = grow_fit_mask(
        coh, vel, seed_threshold=0.85, floor=0.7, similarity=sim, valid=valid, **GROW
    )
    assert not mask[220:226, 180:190].any()
    assert not mask[220:226, 190:200].any()
    assert not mask[220:226, 200:210].any()
    assert mask[220:226, 210:220].all()


def test_grow_excluded_seeds_do_not_set_the_local_velocity():
    """Seeds on a known deforming area stay, but the reference ignores them."""
    coh, vel = _grow_scene()
    coh[180:220, 100:140] = 0.9  # stable seeds at velocity 0
    vel[180:220, 180:220] = -30.0  # subsiding seeds
    coh[220:226, 100:140] = 0.8
    exclude = np.zeros_like(coh, dtype=bool)
    exclude[180:220, 180:220] = True
    mask = grow_fit_mask(
        coh, vel, seed_threshold=0.85, floor=0.7, exclude=exclude, **GROW
    )
    assert mask[180:220, 180:220].all()
    assert mask[220:226, 100:140].all()


def test_grow_rejects_mixed_grids():
    with pytest.raises(ValueError, match="velocity"):
        grow_fit_mask(np.ones((4, 4)), np.ones((4, 5)))


def test_cli_grow_needs_velocity(tmp_path):
    p = _product(tmp_path / "p.nc", np.full((4, 4), 0.9))
    with pytest.raises(SystemExit):
        main([str(p), "--grow-floor", "0.6", "--out", str(tmp_path / "m.tif")])
