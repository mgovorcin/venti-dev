# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Gap filling and residual-outlier masking for raster data."""

from __future__ import annotations

import numpy as np


def fill_gaps(
    array: np.ndarray,
    gaps: np.ndarray,
    max_search_distance: int | None = None,
    smoothing_iterations: int = 0,
) -> np.ndarray:
    """Fill gaps by inverse-distance interpolation from surrounding pixels.

    Uses GDAL's ``FillNodata``. Everything happens in memory: the output
    array is handed to GDAL without copying and GDAL's work files use the
    ``MEM`` driver, so no writable working or temporary directory is needed
    (e.g. in a container). GDAL exceptions are enabled only for this call,
    leaving the process-wide GDAL setting untouched.

    Parameters
    ----------
    array : np.ndarray
        2-D data. Non-finite pixels are treated as gaps too.
    gaps : np.ndarray
        Boolean array, True where pixels should be filled.
    max_search_distance : int, optional
        Maximum distance in pixels to search for valid values, by default
        the largest array dimension. Gaps farther than this stay NaN.
    smoothing_iterations : int, optional
        3x3 smoothing passes over the filled pixels, by default 0.

    Returns
    -------
    np.ndarray
        Filled float32 copy of `array`; valid pixels are unchanged.

    Raises
    ------
    ValueError
        If `gaps` does not match the shape of `array`.

    Examples
    --------
    >>> import numpy as np
    >>> from venti.spatial import fill_gaps
    >>> data = np.array([[1.0, 1.0, 1.0], [1.0, 9.0, 1.0], [1.0, 1.0, 1.0]])
    >>> gaps = data > 5
    >>> fill_gaps(data, gaps)[1, 1]
    np.float32(1.0)

    """
    from osgeo import gdal, gdal_array

    if np.shape(gaps) != np.shape(array):
        msg = f"gaps shape {np.shape(gaps)} != array shape {np.shape(array)}"
        raise ValueError(msg)

    # The only copy: GDAL fills it in place through the array wrapper below.
    filled = np.array(array, dtype=np.float32, order="C")
    valid = ~np.asarray(gaps, dtype=bool) & np.isfinite(filled)
    if valid.all() or not valid.any():
        return filled
    # Gap values are never used as sources; NaN marks any GDAL cannot reach.
    filled[~valid] = np.nan

    if max_search_distance is None:
        max_search_distance = max(filled.shape)
    with gdal.ExceptionMgr(useExceptions=True):
        data_ds = gdal_array.OpenArray(filled)
        mask_ds = gdal_array.OpenArray(valid.view(np.uint8))
        gdal.FillNodata(
            targetBand=data_ds.GetRasterBand(1),
            maskBand=mask_ds.GetRasterBand(1),
            maxSearchDist=max_search_distance,
            smoothingIterations=smoothing_iterations,
            # Work files default to GeoTIFFs in CPL_TMPDIR/TMPDIR or the
            # working directory; in memory they are faster and need no disk.
            options=["TEMP_FILE_DRIVER=MEM"],
        )
        data_ds.FlushCache()
    return filled


def _get_residual_mask(
    insar_data: np.ndarray,
    gnss_los: np.ndarray,
    lower_quantile: float = 0.15,
    upper_quantile: float = 0.85,
) -> np.ndarray:
    """Mask pixels whose InSAR-GNSS residual is an outlier.

    Parameters
    ----------
    insar_data : np.ndarray
        InSAR displacement or velocity field.
    gnss_los : np.ndarray
        GNSS LOS reference field.
    lower_quantile : float, optional
        Lower residual quantile threshold, by default ``0.15``.
    upper_quantile : float, optional
        Upper residual quantile threshold, by default ``0.85``.

    Returns
    -------
    np.ndarray
        Boolean mask: ``True`` where pixels should be excluded.

    """
    invalid = np.ma.masked_invalid(insar_data).mask
    residual_filled = np.where(invalid, np.nan, insar_data - gnss_los)
    lo = np.nanquantile(residual_filled, lower_quantile)
    hi = np.nanquantile(residual_filled, upper_quantile)
    return invalid | (residual_filled < lo) | (residual_filled > hi)
