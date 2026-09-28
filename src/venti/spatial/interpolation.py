"""Spatial interpolation of scattered point data onto regular grids.

Provides RBF and griddata interpolation with tiled evaluation to limit
peak memory usage on large raster grids.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Literal, TypeAlias

import numpy as np
import xarray as xr
from scipy.interpolate import Rbf, RegularGridInterpolator, griddata
from tqdm import tqdm

logger = logging.getLogger(__name__)

InterpolationMethod = Literal["rbf", "griddata"]


def fill_masked_region(
    data: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    """Fill masked pixels with nearest valid neighbor values.

    Pixels where ``mask == 0`` are replaced by the value of their nearest
    spatially valid neighbor — a pixel that is neither masked nor NaN.

    Parameters
    ----------
    data : np.ndarray
        2-D data array.  NaN denotes pixels that are independently invalid
        (e.g. water or no-data regions).
    mask : np.ndarray
        2-D binary mask aligned to ``data`` (1 = valid, 0 = region to fill).

    Returns
    -------
    np.ndarray
        Copy of ``data`` with masked pixels replaced by their nearest valid
        neighbor.

    Notes
    -----
    Uses :func:`scipy.ndimage.distance_transform_edt` to locate the nearest
    valid source pixel for every masked pixel in O(n) time.

    Examples
    --------
    ::

        filled = fill_masked_region(displacement_mm, event_mask)

    """
    from scipy.ndimage import distance_transform_edt

    fill_region = ~mask.astype(bool)
    # Source pixels: not NaN and not in the fill region
    valid = ~np.isnan(data) & ~fill_region

    if not valid.any():
        logger.warning(
            "No valid source pixels outside the masked region; fill has no effect"
        )
        return data.copy()

    # For every pixel in ~valid, find the (row, col) of the nearest valid pixel.
    # With return_distances=False, this returns a single (2, ny, nx) index
    # array directly (not a (distances, indices) tuple) — do not unpack it.
    nearest = distance_transform_edt(
        ~valid, return_distances=False, return_indices=True
    )

    filled = data.copy()
    filled[fill_region] = data[nearest[0][fill_region], nearest[1][fill_region]]
    return filled


def detect_coherent_residual_regions(
    residual: np.ndarray,
    candidate_mad_threshold: float,
    min_region_pixels: int,
) -> np.ndarray:
    """Flag spatially-coherent regions of anomalous residual.

    A global per-pixel MAD test only catches isolated, extreme outliers —
    it structurally misses a wide, gradually-varying feature (e.g. real,
    localized subsidence the GNSS grid's smoothing doesn't resolve), since
    no single pixel in such a region need be extreme relative to the
    frame's overall noise level. This instead requires candidate pixels to
    form a large *connected* region: seed candidates with a loose per-pixel
    threshold, label connected components, and keep only regions at or
    above `min_region_pixels`. Region size, not per-pixel magnitude, is
    what separates real coherent signal from noise here — an isolated bad
    pixel won't form a large enough region even if individually extreme;
    a real deforming area will, even if no single pixel in it crosses a
    strict outlier threshold on its own.

    Parameters
    ----------
    residual : np.ndarray
        2-D residual field (e.g. ``disp - gnss_los``).
    candidate_mad_threshold : float
        Loose per-pixel threshold (in scaled-MAD units, same 1.4826x
        Gaussian consistency correction as a standalone MAD test) used to
        seed candidate pixels before connected-component labeling.
        Deliberately looser than a standalone outlier threshold would be.
    min_region_pixels : int
        Minimum connected-region size (pixels) for a candidate region to
        be flagged as real, coherent deformation.

    Returns
    -------
    np.ndarray
        Boolean mask, ``True`` where pixels belong to a flagged region.

    """
    from scipy import ndimage as ndi

    valid = np.isfinite(residual)
    median_res = np.nanmedian(residual)
    scaled_mad = np.nanmedian(np.abs(residual - median_res)) * 1.4826
    candidate = valid & (
        np.abs(residual - median_res) > candidate_mad_threshold * scaled_mad
    )

    labels, n_labels = ndi.label(candidate)
    if n_labels == 0:
        return np.zeros_like(residual, dtype=bool)

    region_sizes = ndi.sum(candidate, labels, index=np.arange(1, n_labels + 1))
    large_labels = np.where(region_sizes >= min_region_pixels)[0] + 1

    return np.isin(labels, large_labels)


RbfFunction = Literal[
    "multiquadric", "inverse", "gaussian", "linear", "cubic", "quintic", "thin_plate"
]
GriddataMethod = Literal["linear", "nearest", "cubic"]


GridLike: TypeAlias = str | Path | xr.Dataset | tuple[np.ndarray, np.ndarray]


def grid_coordinates(grid: GridLike) -> tuple[np.ndarray, np.ndarray]:
    """Return the 1-D ``(x, y)`` pixel-centre coordinates of a raster grid.

    Parameters
    ----------
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        A NetCDF (``x``/``y`` coordinates, e.g. an OPERA DISP-S1 product), a
        GeoTIFF (from its geotransform), an `xr.Dataset` with ``x``/``y``, or
        an ``(x, y)`` pair of 1-D coordinate arrays.

    Returns
    -------
    x, y : np.ndarray
        Column and row coordinates in the grid's CRS.

    Raises
    ------
    ValueError
        If the grid is rotated or the coordinate arrays are not 1-D.

    Examples
    --------
    >>> import numpy as np
    >>> from venti.spatial.interpolation import grid_coordinates
    >>> x, y = grid_coordinates((np.arange(3.0), np.arange(2.0)))
    >>> x.size, y.size
    (3, 2)

    """
    if isinstance(grid, tuple | list):
        x, y = (np.asarray(v) for v in grid)
    elif isinstance(grid, xr.Dataset):
        x, y = np.asarray(grid["x"].values), np.asarray(grid["y"].values)
    elif Path(grid).suffix.lower() in (".tif", ".tiff"):
        import rasterio

        with rasterio.open(grid) as src:
            t = src.transform
            if t.b != 0 or t.d != 0:
                msg = f"Rotated grids are not supported: {grid}"
                raise ValueError(msg)
            x = t.c + (np.arange(src.width) + 0.5) * t.a
            y = t.f + (np.arange(src.height) + 0.5) * t.e
    else:
        with xr.open_dataset(grid) as ds:
            x, y = np.asarray(ds["x"].values), np.asarray(ds["y"].values)
    if x.ndim != 1 or y.ndim != 1:
        msg = f"Grid coordinates must be 1-D, got x{x.shape}, y{y.shape}"
        raise ValueError(msg)
    return x, y


def _regular_grid_interpolator(
    grid: GridLike, arr: np.ndarray
) -> RegularGridInterpolator:
    """Build a (y, x) regular-grid interpolator for `arr` on `grid`."""
    x, y = grid_coordinates(grid)
    assert arr.shape == (
        y.size,
        x.size,
    ), f"Array shape {arr.shape} must match (len(y)={y.size}, len(x)={x.size})"

    if np.any(np.diff(y) < 0):
        y = y[::-1]
        arr = arr[::-1, :]
    if np.any(np.diff(x) < 0):
        x = x[::-1]
        arr = arr[:, ::-1]

    return RegularGridInterpolator((y, x), arr, bounds_error=False, fill_value=np.nan)


def _sample_on_points(
    grid: GridLike,
    arr: np.ndarray,
    pts_x: np.ndarray,
    pts_y: np.ndarray,
) -> np.ndarray:
    """Sample raster `arr` on `grid` at arbitrary (x, y) points (NaN outside)."""
    interp = _regular_grid_interpolator(grid, arr)
    return interp(np.column_stack([pts_y, pts_x]))


def interpolate_rbf(
    grid: GridLike,
    gx: np.ndarray,
    gy: np.ndarray,
    zv: np.ndarray,
    function: RbfFunction = "cubic",
    tile_nx: int = 512,
    tile_ny: int = 512,
    dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Interpolate scattered samples onto a regular grid using RBF.

    Uses :class:`scipy.interpolate.Rbf` with tiled evaluation to limit
    peak memory usage on large grids.

    Parameters
    ----------
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        Output grid: NetCDF, GeoTIFF, Dataset or ``(x, y)`` coordinate arrays
        (see `grid_coordinates`).
    gx : np.ndarray
        X coordinates of the scattered input samples.
    gy : np.ndarray
        Y coordinates of the scattered input samples.
    zv : np.ndarray
        Values at the scattered input samples.
    function : str, optional
        RBF basis function, by default ``'cubic'``.
    tile_nx : int, optional
        Tile width in pixels, by default ``512``.
    tile_ny : int, optional
        Tile height in pixels, by default ``512``.
    dtype : numpy dtype, optional
        Output dtype, by default ``np.float32``.

    Returns
    -------
    np.ndarray
        Interpolated grid of shape ``(ny, nx)``.

    """
    x, y = (c.astype(dtype) for c in grid_coordinates(grid))

    nx, ny = x.size, y.size
    Z = np.full((ny, nx), np.nan, dtype=dtype)
    rbf = Rbf(gx, gy, zv, function=function, smooth=1)

    nx_tiles = math.ceil(nx / tile_nx)
    ny_tiles = math.ceil(ny / tile_ny)

    with tqdm(total=nx_tiles * ny_tiles, desc="RBF interpolation") as pbar:
        for j in range(ny_tiles):
            y0, y1 = j * tile_ny, min((j + 1) * tile_ny, ny)
            y_slice = y[y0:y1]

            for i in range(nx_tiles):
                x0, x1 = i * tile_nx, min((i + 1) * tile_nx, nx)
                x_slice = x[x0:x1]

                xv, yv = np.meshgrid(x_slice, y_slice, indexing="xy")
                Z[y0:y1, x0:x1] = rbf(xv, yv).reshape(y_slice.size, x_slice.size)
                pbar.update(1)

    return Z


def interpolate_griddata(
    grid: GridLike,
    gx: np.ndarray,
    gy: np.ndarray,
    zv: np.ndarray,
    method: GriddataMethod = "linear",
    tile_nx: int = 512,
    tile_ny: int = 512,
    dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Interpolate scattered samples onto a regular grid using griddata.

    Parameters
    ----------
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        Output grid: NetCDF, GeoTIFF, Dataset or ``(x, y)`` coordinate arrays
        (see `grid_coordinates`).
    gx : np.ndarray
        X coordinates of the scattered input samples.
    gy : np.ndarray
        Y coordinates of the scattered input samples.
    zv : np.ndarray
        Values at the scattered input samples.
    method : str, optional
        Interpolation method passed to :func:`scipy.interpolate.griddata`,
        by default ``'linear'``.
    tile_nx : int, optional
        Tile width in pixels, by default ``512``.
    tile_ny : int, optional
        Tile height in pixels, by default ``512``.
    dtype : numpy dtype, optional
        Output dtype, by default ``np.float32``.

    Returns
    -------
    np.ndarray
        Interpolated grid of shape ``(ny, nx)``.  Unlike RBF, ``'linear'``
        and ``'cubic'`` do not extrapolate beyond the convex hull of the
        input points, which can otherwise leave the output frame's own
        corners/edges as ``NaN`` whenever the input points were themselves
        queried within that same frame's bounds (so the frame's corners sit
        just outside the enclosed points' own hull). Those pixels are
        filled with nearest-neighbor interpolation instead, matching the
        pattern already used for gridded GNSS products in this project's
        demo notebooks — ``method='nearest'`` extrapolates everywhere, so
        this is only ever a fallback for the region the primary `method`
        left uncovered, not a replacement for it.

    """
    x, y = (c.astype(dtype) for c in grid_coordinates(grid))

    nx, ny = x.size, y.size
    Z = np.full((ny, nx), np.nan, dtype=dtype)
    points = np.column_stack((gx, gy))
    values = np.asarray(zv, dtype=dtype)

    nx_tiles = math.ceil(nx / tile_nx)
    ny_tiles = math.ceil(ny / tile_ny)

    with tqdm(total=nx_tiles * ny_tiles, desc="griddata interpolation") as pbar:
        for j in range(ny_tiles):
            y0, y1 = j * tile_ny, min((j + 1) * tile_ny, ny)
            y_slice = y[y0:y1]

            for i in range(nx_tiles):
                x0, x1 = i * tile_nx, min((i + 1) * tile_nx, nx)
                x_slice = x[x0:x1]

                xv, yv = np.meshgrid(x_slice, y_slice, indexing="xy")
                coords = np.column_stack((xv.ravel(), yv.ravel()))
                z_tile = griddata(points, values, coords, method=method)
                Z[y0:y1, x0:x1] = z_tile.reshape(y_slice.size, x_slice.size)
                pbar.update(1)

    if method != "nearest":
        gaps = np.isnan(Z)
        if gaps.any():
            logger.info(
                "griddata (method=%r) left %d px outside the input points' "
                "convex hull; filling with nearest-neighbor interpolation",
                method,
                int(gaps.sum()),
            )
            xv, yv = np.meshgrid(x, y, indexing="xy")
            gap_coords = np.column_stack((xv[gaps], yv[gaps]))
            Z[gaps] = griddata(points, values, gap_coords, method="nearest")

    return Z
