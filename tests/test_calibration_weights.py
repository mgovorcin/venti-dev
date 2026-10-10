# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Fit weights: coherence^p x local robust weights (PRD R-S4, plan T31)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.calibration.loclin import kernel_sigma_px, local_linear_surface
from venti.calibration.weights import coherence_weights, fit_weights, robust_weights

PIXEL_M = 180.0
SIGMA = kernel_sigma_px(50_000.0, PIXEL_M)


def test_coherence_weights():
    coh = np.array([[0.0, 0.5, 1.0, np.nan], [1.2, 0.9, -0.1, 0.25]])
    w = coherence_weights(coh, 8.0, coh.shape)
    assert w[0, 0] == 0.0
    assert w[0, 1] == pytest.approx(0.5**8)
    assert w[0, 2] == 1.0
    assert w[0, 3] == 0.0  # NaN -> 0
    assert w[1, 0] == 1.0  # clipped
    assert w[1, 2] == 0.0
    np.testing.assert_array_equal(coherence_weights(coh, 0.0, coh.shape), 1.0)
    np.testing.assert_array_equal(coherence_weights(None, 8.0, (2, 4)), 1.0)
    with pytest.raises(ValueError, match="does not match"):
        coherence_weights(coh, 2.0, (3, 3))
    with pytest.raises(ValueError, match=">= 0"):
        coherence_weights(coh, -1.0, coh.shape)


def _field_with_outlier(outlier_mm: float):
    rng = np.random.default_rng(0)
    ny, nx = 160, 200
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    truth = 3.0 + 0.01 * xx - 0.005 * yy
    field = truth + 0.5 * rng.standard_normal((ny, nx))
    field[78:82, 98:102] += outlier_mm  # a 4x4 px blunder
    return field, truth


@pytest.mark.parametrize("method", ["huber", "gate"])
def test_robust_weights_flag_the_outlier_only(method):
    field, _ = _field_with_outlier(100.0)
    base = np.ones_like(field)
    rw = robust_weights(field, base, SIGMA, method=method)
    assert rw[78:82, 98:102].max() == 0.0  # 100 mm vs 0.5 mm noise: gated out
    clean = np.ones_like(rw, bool)
    clean[70:90, 90:110] = False
    if method == "gate":
        assert np.mean(rw[clean] == 1.0) > 0.98  # honest pixels keep weight 1
    else:
        # Huber trims the ~18 % of Gaussian residuals beyond 1.345 MAD a little
        assert np.mean(rw[clean]) > 0.9
        assert np.mean(rw[clean] > 0.5) > 0.98
    assert rw.min() >= 0
    assert rw.max() <= 1


def test_single_outlier_does_not_move_the_surface():
    """A 10 cm blunder must move the 50 km surface by < 0.1 mm once weighted (R-S4)."""
    clean_field, truth = _field_with_outlier(0.0)
    dirty_field, _ = _field_with_outlier(100.0)
    base = np.ones_like(truth)
    s_clean, _ = local_linear_surface(clean_field, base, SIGMA)
    w = fit_weights(dirty_field, base, SIGMA, robust=True)
    s_robust, _ = local_linear_surface(dirty_field, w, SIGMA)
    s_naive, _ = local_linear_surface(dirty_field, base, SIGMA)
    assert np.abs(s_robust - s_clean).max() < 0.1
    assert np.abs(s_naive - s_clean).max() > 0.1  # the blunder does matter unweighted


def test_fit_weights_defaults_are_the_identity():
    field = np.random.default_rng(1).standard_normal((30, 40))
    base = np.random.default_rng(2).uniform(0, 1, (30, 40))
    np.testing.assert_array_equal(fit_weights(field, base, 5.0), base)


def test_fit_weights_combine_coherence_and_robust():
    field, _ = _field_with_outlier(100.0)
    base = np.ones_like(field)
    coh = np.full_like(field, 0.9)
    coh[:, :100] = 0.6
    w_coh = fit_weights(field, base, SIGMA, coherence=coh, coherence_power=8)
    assert w_coh[10, 10] == pytest.approx(0.6**8)
    assert w_coh[10, 150] == pytest.approx(0.9**8)
    w = fit_weights(field, base, SIGMA, coherence=coh, coherence_power=8, robust=True)
    assert w[78:82, 98:102].max() == 0.0
    assert np.all(w <= w_coh + 1e-12)  # robust only ever lowers a weight


def test_robust_weights_respect_zero_base_weight_and_validate():
    field, _ = _field_with_outlier(0.0)
    base = np.ones_like(field)
    base[:, :20] = 0.0
    rw = robust_weights(field, base, SIGMA)
    assert rw[:, :20].max() == 0.0
    with pytest.raises(ValueError, match="differ"):
        robust_weights(field, base[:10], SIGMA)
    with pytest.raises(ValueError, match="iterations"):
        robust_weights(field, base, SIGMA, iterations=0)
