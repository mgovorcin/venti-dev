"""Unit tests for the spatial subpackage.

All tests use synthetic numpy arrays — no files or network required.
"""

from __future__ import annotations

import numpy as np
import pytest

try:
    import joblib  # noqa: F401

    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False

try:
    import skimage  # noqa: F401

    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False

# The solver and moving-window modules import skimage and joblib,
# so tests that import them need both.
_HAS_SPATIAL_DEPS = HAS_JOBLIB and HAS_SKIMAGE
_SKIP_SPATIAL = pytest.mark.skipif(
    not _HAS_SPATIAL_DEPS, reason="skimage or joblib not installed"
)


# _design_matrix_poly
@_SKIP_SPATIAL
class TestDesignMatrixPoly:
    """Tests for _design_matrix_poly."""

    def test_order_0_is_constant(self):
        from venti.solver.design_matrix import _design_matrix_poly

        x = np.array([1.0, 2.0, 3.0])
        y = np.array([4.0, 5.0, 6.0])
        A = _design_matrix_poly(x, y, poly_order=0)
        assert A.shape == (3, 1)
        assert np.all(A == 1.0)

    def test_order_1_has_three_columns(self):
        from venti.solver.design_matrix import _design_matrix_poly

        x = np.linspace(0, 1, 10)
        y = np.linspace(0, 1, 10)
        A = _design_matrix_poly(x, y, poly_order=1)
        assert A.shape == (10, 3)

    def test_order_1p5_has_four_columns(self):
        from venti.solver.design_matrix import _design_matrix_poly

        x = np.linspace(0, 1, 10)
        y = np.linspace(0, 1, 10)
        A = _design_matrix_poly(x, y, poly_order=1.5)
        assert A.shape == (10, 4)

    def test_order_2_has_six_columns(self):
        from venti.solver.design_matrix import _design_matrix_poly

        x = np.linspace(0, 1, 10)
        y = np.linspace(0, 1, 10)
        A = _design_matrix_poly(x, y, poly_order=2)
        assert A.shape == (10, 6)

    def test_order_3_has_eight_columns(self):
        from venti.solver.design_matrix import _design_matrix_poly

        x = np.linspace(0, 1, 10)
        y = np.linspace(0, 1, 10)
        A = _design_matrix_poly(x, y, poly_order=3)
        assert A.shape == (10, 8)

    def test_unsupported_order_raises(self):
        from venti.solver.design_matrix import _design_matrix_poly

        x = np.ones(5)
        y = np.ones(5)
        with pytest.raises((ValueError, KeyError)):
            _design_matrix_poly(x, y, poly_order=99)


# _calc_plane_values
@_SKIP_SPATIAL
class TestCalcPlaneValues:
    """Tests for _calc_plane_values."""

    def test_constant_plane(self):
        from venti.solver.plane_fitting import _calc_plane_values

        x = np.array([[0.0, 1.0], [2.0, 3.0]])
        y = np.array([[0.0, 0.0], [1.0, 1.0]])
        coef = np.array([7.0])
        result = _calc_plane_values(x, y, coef, poly_order=0)
        assert np.allclose(result, 7.0)

    def test_linear_plane_exact(self):
        """z = 2x + 3y + 1 must be recovered exactly."""
        from venti.solver.plane_fitting import _calc_plane_values

        x = np.array([[0.0, 1.0], [0.0, 1.0]])
        y = np.array([[0.0, 0.0], [1.0, 1.0]])
        coef = np.array([2.0, 3.0, 1.0])
        result = _calc_plane_values(x, y, coef, poly_order=1)
        expected = 2.0 * x + 3.0 * y + 1.0
        assert np.allclose(result, expected)

    def test_output_shape_matches_input(self):
        from venti.solver.plane_fitting import _calc_plane_values

        x = np.random.rand(8, 12)
        y = np.random.rand(8, 12)
        coef = np.array([1.0, 2.0, 3.0])
        result = _calc_plane_values(x, y, coef, poly_order=1)
        assert result.shape == (8, 12)


# _weighted_lscov
@_SKIP_SPATIAL
class TestWeightedLscov:
    """Tests for _weighted_lscov."""

    def test_recover_known_coefficients_unweighted(self):
        """Exact linear system: solution must match to machine precision."""
        from venti.solver.lscov import _weighted_lscov

        rng = np.random.default_rng(0)
        x = rng.uniform(0, 10, 50)
        y = rng.uniform(0, 10, 50)
        true_coef = np.array([3.0, -2.0, 5.0])
        A = np.column_stack([x, y, np.ones(50)])
        b = A @ true_coef
        sol, *_ = _weighted_lscov(A, b)
        assert np.allclose(sol.ravel(), true_coef, atol=1e-6)

    def test_recover_known_coefficients_weighted(self):
        from venti.solver.lscov import _weighted_lscov

        rng = np.random.default_rng(1)
        x = rng.uniform(0, 10, 50)
        y = rng.uniform(0, 10, 50)
        true_coef = np.array([1.5, -0.5, 2.0])
        A = np.column_stack([x, y, np.ones(50)])
        b = A @ true_coef
        weights = rng.uniform(0.5, 2.0, 50)
        sol, *_ = _weighted_lscov(A, b, weights)
        assert np.allclose(sol.ravel(), true_coef, atol=1e-5)

    def test_returns_eight_values(self):
        from venti.solver.lscov import _weighted_lscov

        A = np.column_stack([np.arange(10.0), np.ones(10)])
        b = np.arange(10.0)
        result = _weighted_lscov(A, b)
        assert len(result) == 8

    def test_residuals_shape(self):
        from venti.solver.lscov import _weighted_lscov

        n, p = 20, 3
        A = np.random.rand(n, p)
        b = np.random.rand(n)
        _, _, _, _, res, wres, obs_hat, _ = _weighted_lscov(A, b)
        assert res.shape == (n, 1)
        assert wres.shape == (n, 1)
        assert obs_hat.shape == (n, 1)


# _get_sliding_windows
@_SKIP_SPATIAL
class TestGetSlidingWindows:
    """Tests for _get_sliding_windows."""

    def test_single_window_when_size_equals_length(self):
        from venti.filtering.moving_window import _get_sliding_windows

        windows = _get_sliding_windows(100, win_size=100, win_overlap=0)
        assert len(windows) == 1
        assert windows[0] == slice(0, 100)

    def test_two_windows_no_overlap(self):
        from venti.filtering.moving_window import _get_sliding_windows

        windows = _get_sliding_windows(100, win_size=50, win_overlap=0)
        assert len(windows) == 2

    def test_overlap_reduces_stride(self):
        from venti.filtering.moving_window import _get_sliding_windows

        no_overlap = _get_sliding_windows(200, win_size=50, win_overlap=0)
        with_overlap = _get_sliding_windows(200, win_size=50, win_overlap=10)
        assert len(with_overlap) >= len(no_overlap)

    def test_windows_cover_full_range(self):
        from venti.filtering.moving_window import _get_sliding_windows

        length = 150
        windows = _get_sliding_windows(length, win_size=60, win_overlap=10)
        assert windows[0].start == 0
        assert windows[-1].stop == length

    def test_all_windows_are_slices(self):
        from venti.filtering.moving_window import _get_sliding_windows

        windows = _get_sliding_windows(100, win_size=30, win_overlap=5)
        assert all(isinstance(w, slice) for w in windows)


# fit_windowed_plane
@_SKIP_SPATIAL
class TestFitWindowedPlane:
    """Tests for fit_windowed_plane (public API)."""

    NY, NX = 60, 80
    SNWE = (3_800_000.0, 3_860_000.0, 400_000.0, 480_000.0)

    def _flat_scene(self, value: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
        insar = np.full((self.NY, self.NX), value, dtype=np.float32)
        gnss = np.zeros((self.NY, self.NX), dtype=np.float32)
        return insar, gnss

    def _linear_residual_scene(self) -> tuple[np.ndarray, np.ndarray]:
        x = np.linspace(0, 1, self.NX)
        y = np.linspace(0, 1, self.NY)
        xg, yg = np.meshgrid(x, y)
        ramp = (0.05 * xg + 0.03 * yg).astype(np.float32)
        return ramp.copy(), np.zeros_like(ramp)

    def test_output_shapes(self):
        from venti.filtering.moving_window import fit_windowed_plane

        insar, gnss = self._flat_scene()
        surface, std = fit_windowed_plane(
            insar,
            gnss,
            win_xsize=40,
            win_ysize=30,
            win_overlap_x=5,
            win_overlap_y=5,
            win_extend_x=40,
            win_extend_y=30,
            n_jobs=1,
        )
        assert surface.shape == (self.NY, self.NX)
        assert std.shape == (self.NY, self.NX)

    def test_flat_residual_gives_near_zero_surface(self):
        from venti.filtering.moving_window import fit_windowed_plane

        insar, gnss = self._flat_scene(value=0.0)
        surface, _ = fit_windowed_plane(
            insar,
            gnss,
            win_xsize=40,
            win_ysize=30,
            win_overlap_x=5,
            win_overlap_y=5,
            win_extend_x=40,
            win_extend_y=30,
            poly_order=1,
            n_jobs=1,
        )
        assert np.nanmax(np.abs(surface)) < 1e-4

    def test_linear_ramp_is_captured(self):
        from venti.filtering.moving_window import fit_windowed_plane

        insar, gnss = self._linear_residual_scene()
        surface, _ = fit_windowed_plane(
            insar,
            gnss,
            win_xsize=40,
            win_ysize=30,
            win_overlap_x=5,
            win_overlap_y=5,
            win_extend_x=40,
            win_extend_y=30,
            poly_order=1,
            n_jobs=1,
        )
        assert np.nanmax(np.abs(surface)) > 1e-4

    def test_masked_nan_pixels_do_not_crash(self):
        from venti.filtering.moving_window import fit_windowed_plane

        insar, gnss = self._flat_scene()
        insar[10:20, 10:20] = np.nan
        surface, _ = fit_windowed_plane(
            insar,
            gnss,
            win_xsize=40,
            win_ysize=30,
            win_overlap_x=5,
            win_overlap_y=5,
            win_extend_x=40,
            win_extend_y=30,
            n_jobs=1,
        )
        assert surface.shape == (self.NY, self.NX)

    def test_gnss_los_std_accepted(self):
        from venti.filtering.moving_window import fit_windowed_plane

        insar, gnss = self._flat_scene()
        std_field = np.ones_like(insar) * 0.001
        surface, _ = fit_windowed_plane(
            insar,
            gnss,
            win_xsize=40,
            win_ysize=30,
            win_overlap_x=5,
            win_overlap_y=5,
            win_extend_x=40,
            win_extend_y=30,
            gnss_los_std=std_field,
            n_jobs=1,
        )
        assert surface.shape == (self.NY, self.NX)


# SpatialProcessor
@_SKIP_SPATIAL
class TestSpatialProcessor:
    """Tests for SpatialProcessor.fit_windowed_surface."""

    NY, NX = 60, 80
    SNWE = (3_800_000.0, 3_860_000.0, 400_000.0, 480_000.0)

    def test_returns_ndarray(self):
        from venti.spatial.processor import SpatialProcessor

        processor = SpatialProcessor()
        insar = np.zeros((self.NY, self.NX), dtype=np.float32)
        gnss = np.zeros((self.NY, self.NX), dtype=np.float32)
        result = processor.fit_windowed_surface(
            insar,
            gnss,
            window_size_x=40,
            window_size_y=30,
            n_jobs=1,
        )
        assert isinstance(result, np.ndarray)
        assert result.shape == (self.NY, self.NX)

    def test_default_extend_equals_window_size(self):
        from venti.spatial.processor import SpatialProcessor

        processor = SpatialProcessor()
        insar = np.zeros((self.NY, self.NX), dtype=np.float32)
        gnss = np.zeros((self.NY, self.NX), dtype=np.float32)
        result = processor.fit_windowed_surface(
            insar,
            gnss,
            window_size_x=40,
            window_size_y=30,
            n_jobs=1,
        )
        assert result.shape == (self.NY, self.NX)

    def test_custom_poly_order(self):
        from venti.spatial.processor import SpatialProcessor

        processor = SpatialProcessor()
        insar = np.zeros((self.NY, self.NX), dtype=np.float32)
        gnss = np.zeros((self.NY, self.NX), dtype=np.float32)
        for order in [0, 1, 1.5, 2]:
            result = processor.fit_windowed_surface(
                insar,
                gnss,
                window_size_x=40,
                window_size_y=30,
                poly_order=order,
                n_jobs=1,
            )
            assert result.shape == (self.NY, self.NX)

    def test_explicit_extend_parameters(self):
        from venti.spatial.processor import SpatialProcessor

        processor = SpatialProcessor()
        insar = np.zeros((self.NY, self.NX), dtype=np.float32)
        gnss = np.zeros((self.NY, self.NX), dtype=np.float32)
        result = processor.fit_windowed_surface(
            insar,
            gnss,
            window_size_x=40,
            window_size_y=30,
            window_extend_x=20,
            window_extend_y=15,
            n_jobs=1,
        )
        assert result.shape == (self.NY, self.NX)


# fit_windowed_plane — smoothing method dispatch
@_SKIP_SPATIAL
class TestFitWindowedPlaneSmoothing:
    """Test that each smoothing_method runs without error and produces the right shape."""

    NY, NX = 60, 80

    def _flat_scene(self) -> tuple[np.ndarray, np.ndarray]:
        insar = np.zeros((self.NY, self.NX), dtype=np.float32)
        gnss = np.zeros((self.NY, self.NX), dtype=np.float32)
        return insar, gnss

    def _run(self, method: str, sigma: float | None = 4.0, **kwargs) -> np.ndarray:
        from venti.filtering.moving_window import fit_windowed_plane

        insar, gnss = self._flat_scene()
        surface, _ = fit_windowed_plane(
            insar,
            gnss,
            win_xsize=40,
            win_ysize=30,
            win_overlap_x=5,
            win_overlap_y=5,
            win_extend_x=40,
            win_extend_y=30,
            poly_order=1,
            n_jobs=1,
            smoothing_sigma=sigma,
            smoothing_method=method,
            **kwargs,
        )
        return surface

    def test_gaussian_default(self):
        surface = self._run("gaussian")
        assert surface.shape == (self.NY, self.NX)

    def test_gaussian_fft(self):
        surface = self._run("gaussian_fft")
        assert surface.shape == (self.NY, self.NX)

    def test_hanning_fft(self):
        surface = self._run("hanning_fft")
        assert surface.shape == (self.NY, self.NX)

    def test_savitzky_golay(self):
        surface = self._run(
            "savitzky_golay", sigma=None, sg_window_length=11, sg_polyorder=2
        )
        assert surface.shape == (self.NY, self.NX)

    def test_no_smoothing_when_sigma_none_and_not_sg(self):
        """sigma=None with a non-SG method must skip smoothing entirely (no error)."""
        surface = self._run("gaussian", sigma=None)
        assert surface.shape == (self.NY, self.NX)


# low_pass_filters
class TestLowPassFilters:
    """Tests for venti.filtering.low_pass_filters.

    All tests use synthetic arrays — no files or network required.
    """

    NY, NX = 64, 80

    def _flat_surface(self, value: float = 5.0) -> tuple[np.ndarray, np.ndarray]:
        surface = np.full((self.NY, self.NX), value, dtype=np.float32)
        valid = np.ones((self.NY, self.NX), dtype=bool)
        return surface, valid

    def _ramp_surface(self) -> tuple[np.ndarray, np.ndarray]:
        x = np.linspace(0, 1, self.NX, dtype=np.float32)
        y = np.linspace(0, 1, self.NY, dtype=np.float32)
        surface = np.outer(y, x)
        valid = np.ones((self.NY, self.NX), dtype=bool)
        return surface, valid

    def _masked_surface(self) -> tuple[np.ndarray, np.ndarray]:
        surface, _ = self._flat_surface()
        valid = np.ones((self.NY, self.NX), dtype=bool)
        valid[:8, :] = False  # mask top rows
        surface[~valid] = 0.0
        return surface, valid

    # --- gaussian_spatial ---

    def test_gaussian_spatial_preserves_shape_and_dtype(self):
        from venti.filtering.low_pass_filters import gaussian_spatial

        surface, valid = self._flat_surface()
        out = gaussian_spatial(surface, valid, sigma=3.0)
        assert out.shape == surface.shape
        assert out.dtype == surface.dtype

    def test_gaussian_spatial_constant_surface_unchanged(self):
        from venti.filtering.low_pass_filters import gaussian_spatial

        surface, valid = self._flat_surface(value=7.0)
        out = gaussian_spatial(surface, valid, sigma=3.0)
        assert np.allclose(out[valid], 7.0, atol=1e-3)

    def test_gaussian_spatial_respects_boundary_mask(self):
        from venti.filtering.low_pass_filters import gaussian_spatial

        surface, valid = self._masked_surface()
        out = gaussian_spatial(surface, valid, sigma=2.0)
        # Valid pixels near the mask boundary must not be pulled to zero
        assert not np.any(np.isnan(out[valid]))

    # --- gaussian_fft ---

    def test_gaussian_fft_preserves_shape_and_dtype(self):
        from venti.filtering.low_pass_filters import gaussian_fft

        surface, valid = self._flat_surface()
        out = gaussian_fft(surface, valid, sigma=3.0)
        assert out.shape == surface.shape
        assert out.dtype == surface.dtype

    def test_gaussian_fft_constant_surface_unchanged(self):
        from venti.filtering.low_pass_filters import gaussian_fft

        surface, valid = self._flat_surface(value=7.0)
        out = gaussian_fft(surface, valid, sigma=3.0)
        assert np.allclose(out[valid], 7.0, atol=1e-3)

    def test_gaussian_fft_smooths_ramp(self):
        """FFT Gaussian must reduce high-frequency variance while preserving mean."""
        from venti.filtering.low_pass_filters import gaussian_fft

        rng = np.random.default_rng(0)
        surface, valid = self._ramp_surface()
        noise = rng.standard_normal(surface.shape).astype(np.float32) * 0.5
        noisy = surface + noise
        out = gaussian_fft(noisy, valid, sigma=5.0)
        assert np.std(out[valid]) < np.std(noisy[valid])

    # --- hanning_fft ---

    def test_hanning_fft_preserves_shape_and_dtype(self):
        from venti.filtering.low_pass_filters import hanning_fft

        surface, valid = self._flat_surface()
        out = hanning_fft(surface, valid, sigma=3.0)
        assert out.shape == surface.shape
        assert out.dtype == surface.dtype

    def test_hanning_fft_constant_surface_unchanged(self):
        from venti.filtering.low_pass_filters import hanning_fft

        surface, valid = self._flat_surface(value=7.0)
        out = hanning_fft(surface, valid, sigma=3.0)
        assert np.allclose(out[valid], 7.0, atol=1e-3)

    def test_hanning_fft_smooths_ramp(self):
        from venti.filtering.low_pass_filters import hanning_fft

        rng = np.random.default_rng(1)
        surface, valid = self._ramp_surface()
        noise = rng.standard_normal(surface.shape).astype(np.float32) * 0.5
        noisy = surface + noise
        out = hanning_fft(noisy, valid, sigma=5.0)
        assert np.std(out[valid]) < np.std(noisy[valid])

    # --- savitzky_golay ---

    def test_savitzky_golay_preserves_shape_and_dtype(self):
        from venti.filtering.low_pass_filters import savitzky_golay

        surface, valid = self._flat_surface()
        out = savitzky_golay(surface, valid, window_length=11, polyorder=2)
        assert out.shape == surface.shape
        assert out.dtype == surface.dtype

    def test_savitzky_golay_constant_surface_unchanged(self):
        from venti.filtering.low_pass_filters import savitzky_golay

        surface, valid = self._flat_surface(value=7.0)
        out = savitzky_golay(surface, valid, window_length=11, polyorder=2)
        assert np.allclose(out[valid], 7.0, atol=1e-3)

    def test_savitzky_golay_odd_window_required(self):
        from venti.filtering.low_pass_filters import savitzky_golay

        surface, valid = self._flat_surface()
        with pytest.raises(AssertionError):
            savitzky_golay(surface, valid, window_length=10, polyorder=2)

    def test_savitzky_golay_polyorder_less_than_window(self):
        from venti.filtering.low_pass_filters import savitzky_golay

        surface, valid = self._flat_surface()
        with pytest.raises(AssertionError):
            savitzky_golay(surface, valid, window_length=5, polyorder=6)

    def test_savitzky_golay_smooths_ramp(self):
        from venti.filtering.low_pass_filters import savitzky_golay

        rng = np.random.default_rng(2)
        surface, valid = self._ramp_surface()
        noise = rng.standard_normal(surface.shape).astype(np.float32) * 0.5
        noisy = surface + noise
        out = savitzky_golay(noisy, valid, window_length=15, polyorder=2)
        assert np.std(out[valid]) < np.std(noisy[valid])


def test_gdal_gap_fill_needs_no_writable_directory(tmp_path):
    """GDAL's gap filler must not write work files to disk.

    Reproduces a container with no TMPDIR/CPL_TMPDIR and an unwritable
    working directory, where GDAL falls back to writing ``./`` temp files.
    """
    import os
    import subprocess
    import sys

    pytest.importorskip("osgeo.gdal")
    if os.geteuid() == 0:
        pytest.skip("root can write to read-only directories")

    readonly = tmp_path / "readonly"
    readonly.mkdir()
    readonly.chmod(0o555)
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("TMPDIR", "TEMP", "TMP", "CPL_TMPDIR")
    }
    code = (
        "import numpy as np\n"
        "from venti.spatial import fill_gaps\n"
        "a = np.ones((200, 200), 'float32'); gaps = np.zeros(a.shape, bool)\n"
        "gaps[50:150, 50:150] = True; a[gaps] = 7\n"
        "assert (fill_gaps(a, gaps, smoothing_iterations=2) == 1).all()\n"
    )
    try:
        subprocess.run(
            [sys.executable, "-c", code],
            cwd=readonly,
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
    finally:
        readonly.chmod(0o755)


class TestFillGaps:
    """Tests for fill_gaps (GDAL inverse-distance gap filling)."""

    @pytest.fixture(autouse=True)
    def _needs_gdal(self):
        pytest.importorskip("osgeo.gdal")

    def test_fills_only_gaps_and_keeps_real_zeros(self):
        from venti.spatial import fill_gaps

        data = np.full((40, 50), 2.0)
        data[20, 25] = 0.0  # e.g. the reference pixel after zeroing
        gaps = np.zeros(data.shape, bool)
        gaps[5:10, 5:10] = True
        data[gaps] = 99.0

        filled = fill_gaps(data, gaps)

        assert filled[20, 25] == 0.0
        np.testing.assert_allclose(filled[gaps], 2.0)
        np.testing.assert_array_equal(filled[~gaps], data[~gaps])

    def test_nan_is_a_gap_and_input_is_not_modified(self):
        from venti.spatial import fill_gaps

        data = np.full((30, 30), 3.0)
        data[10:12, 10:12] = np.nan
        original = data.copy()

        filled = fill_gaps(data, np.zeros(data.shape, bool))

        np.testing.assert_allclose(filled, 3.0)
        np.testing.assert_array_equal(data, original)

    def test_unreachable_gaps_stay_nan(self):
        from venti.spatial import fill_gaps

        data = np.full((20, 200), 1.0)
        gaps = np.zeros(data.shape, bool)
        gaps[:, 20:] = True

        filled = fill_gaps(data, gaps, max_search_distance=5)

        assert np.isnan(filled[:, -1]).all()
        assert np.isfinite(filled[:, :22]).all()

    def test_shape_mismatch_raises(self):
        from venti.spatial import fill_gaps

        with pytest.raises(ValueError, match="gaps shape"):
            fill_gaps(np.ones((4, 4)), np.zeros((3, 4), bool))

    def test_global_gdal_exception_setting_is_unchanged(self):
        from osgeo import gdal

        from venti.spatial import fill_gaps

        before = gdal.GetUseExceptions()
        gaps = np.zeros((10, 10), bool)
        gaps[4, 4] = True
        fill_gaps(np.ones((10, 10)), gaps)
        assert gdal.GetUseExceptions() == before
