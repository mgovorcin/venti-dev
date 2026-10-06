# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Gap filling before the fit (PRD R-S3, plan T29)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.calibration import base_weights, fill_gaps

METHODS = ["nearest", "nearest_smooth", "biharmonic", "idw"]


def holes(shape=(40, 50), seed=0, frac=0.3, block=True):
    rng = np.random.default_rng(seed)
    mask = rng.random(shape) > frac
    if block:
        mask[10:20, 15:30] = False  # an enclosed gap
        mask[:, :3] = False  # an edge gap (water along one side)
    return mask


@pytest.mark.parametrize("method", METHODS)
def test_constant_field_is_recovered_exactly(method):
    data = np.full((40, 50), 7.5)
    mask = holes()
    data[~mask] = np.nan
    filled, filled_mask = fill_gaps(data, mask, method=method)
    assert np.isfinite(filled).all()
    np.testing.assert_allclose(filled, 7.5, atol=1e-9)
    np.testing.assert_array_equal(filled_mask, ~mask)


@pytest.mark.parametrize("method", METHODS)
def test_ramp_is_recovered_within_tolerance(method):
    yy, xx = np.mgrid[0:40, 0:50].astype(float)
    truth = 0.4 * xx - 0.2 * yy + 3.0
    mask = holes()
    data = np.where(mask, truth, np.nan)
    filled, _ = fill_gaps(data, mask, method=method)
    # valid pixels untouched
    np.testing.assert_array_equal(filled[mask], truth[mask])
    err = np.abs(filled - truth)[~mask]
    # biharmonic inpainting is for enclosed gaps; it extrapolates poorly into the
    # edge gap, hence its looser bound (the enclosed gap is tested separately)
    tol = {"nearest": 3.0, "nearest_smooth": 2.5, "biharmonic": 2.5, "idw": 1.5}[method]
    if method == "biharmonic":
        enclosed = np.zeros_like(mask)
        enclosed[10:20, 15:30] = True
        assert np.abs(filled - truth)[enclosed].max() < 0.6
    assert err.max() < tol, (method, err.max())


def test_valid_pixels_are_never_changed_and_no_zeros_appear():
    rng = np.random.default_rng(3)
    data = 5.0 + rng.standard_normal((30, 30))
    mask = holes((30, 30), seed=5, frac=0.6)
    data[~mask] = np.nan
    filled, filled_mask = fill_gaps(data, mask)
    np.testing.assert_array_equal(filled[mask], data[mask])
    assert np.isfinite(filled).all()
    assert not np.any(filled == 0.0)  # the gamma failure mode
    assert filled_mask.sum() == (~mask).sum()


def test_nan_in_data_counts_as_invalid_even_if_mask_says_valid():
    data = np.ones((10, 10))
    data[4, 4] = np.nan
    filled, filled_mask = fill_gaps(data, np.ones((10, 10), bool))
    assert filled[4, 4] == 1.0
    assert filled_mask[4, 4]
    assert filled_mask.sum() == 1


def test_nothing_to_fill_returns_copy():
    data = np.arange(12.0).reshape(3, 4)
    filled, filled_mask = fill_gaps(data)
    np.testing.assert_array_equal(filled, data)
    assert filled is not data
    assert not filled_mask.any()


def test_all_invalid_raises():
    with pytest.raises(ValueError, match="no valid pixel"):
        fill_gaps(np.full((5, 5), np.nan))
    with pytest.raises(ValueError, match="2-D"):
        fill_gaps(np.ones(5))


@pytest.mark.parametrize("seed", range(5))
def test_property_random_masks_give_finite_fields(seed):
    rng = np.random.default_rng(seed)
    shape = (25, 33)
    data = rng.standard_normal(shape)
    mask = rng.random(shape) > rng.uniform(0.2, 0.95)
    mask.flat[rng.integers(0, mask.size)] = True  # at least one valid pixel
    data[~mask] = np.nan
    filled, filled_mask = fill_gaps(data, mask)
    assert np.isfinite(filled).all()
    assert filled_mask.sum() == (~mask).sum()


def test_base_weights():
    valid = np.array([[True, False, False], [True, True, False]])
    filled = np.array([[False, True, False], [False, False, True]])
    w = base_weights(valid, filled, w_filled=0.02)
    np.testing.assert_array_equal(w, [[1.0, 0.02, 0.0], [1.0, 1.0, 0.02]])
    # filled wins nothing over valid
    w2 = base_weights(valid, valid | filled)
    assert w2[0, 0] == 1.0
    with pytest.raises(ValueError, match="differ"):
        base_weights(valid, filled[:1])
    with pytest.raises(ValueError, match="w_filled"):
        base_weights(valid, filled, w_filled=2.0)


def test_weights_from_fill_gaps_output():
    data = np.ones((8, 8))
    mask = np.ones((8, 8), bool)
    mask[2:4, 2:4] = False
    data[~mask] = np.nan
    _filled, filled_mask = fill_gaps(data, mask)
    w = base_weights(mask, filled_mask)
    assert w.sum() == pytest.approx(60 + 4 * 0.02)
