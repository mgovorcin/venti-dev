# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""LOS projection for GNSS displacement data.

Projects GNSS east/north/up observations into the InSAR line-of-sight
direction and interpolates the result onto the full raster grid.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from ..spatial.interpolation import (
    GridLike,
    InterpolationMethod,
    RbfFunction,
    _sample_on_points,
    grid_coordinates,
    interpolate_griddata,
    interpolate_rbf,
)

if TYPE_CHECKING:
    import geopandas as gpd

logger = logging.getLogger(__name__)


def _sample_los(coords, los_east, los_north, los_up, x, y, dtype):
    """Sample the LOS unit vectors at points; NaN where there is no look.

    LOS rasters mark pixels outside the swath either with NaN or with
    ``(e, n, u) = (0, 0, 1)`` / ``(0, 0, 0)``. A SAR look is never
    vertical, so ``e == n == 0`` is no look at all; using it would enter
    such a station into the GNSS field with its vertical component only.
    """
    e, n, u = (
        _sample_on_points(coords, np.asarray(c, dtype=dtype), x, y)
        for c in (los_east, los_north, los_up)
    )
    no_look = (e == 0) & (n == 0)
    return tuple(np.where(no_look, np.nan, c) for c in (e, n, u))


def project_to_los(
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    grid: GridLike,
    gnss_gdf: gpd.GeoDataFrame,
    method: InterpolationMethod = "rbf",
    rbf_function: RbfFunction = "cubic",
) -> np.ndarray:
    """Project GNSS east/north/up displacements into the InSAR LOS direction.

    Samples the LOS unit vectors at each GNSS station location, forms the
    dot product with the GNSS displacement vector, then interpolates the
    resulting scalar field onto the full raster grid. Stations where the
    LOS is NaN or has no horizontal component (the ``(0, 0, 1)`` fill some
    LOS rasters use outside the swath) are left out.

    Parameters
    ----------
    los_east : np.ndarray
        2-D array of east LOS unit-vector components, shape ``(ny, nx)``.
    los_north : np.ndarray
        2-D array of north LOS unit-vector components.
    los_up : np.ndarray
        2-D array of up LOS unit-vector components.
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        Grid of the LOS arrays and the output: NetCDF, GeoTIFF, Dataset or
        ``(x, y)`` coordinate arrays, in the stations' CRS (see
        `venti.spatial.interpolation.grid_coordinates`).
    gnss_gdf : gpd.GeoDataFrame
        GNSS station data with a ``geometry`` column (Point, UTM) and columns
        ``'deast'``, ``'dnorth'``, ``'dup'`` containing the displacement or
        velocity in each component.
    method : {'rbf', 'griddata'}, optional
        Spatial interpolation method, by default ``'rbf'``.
    rbf_function : str, optional
        RBF kernel passed to :func:`~venti.spatial.interpolation.interpolate_rbf`
        when ``method='rbf'``, by default ``'cubic'``.

    Returns
    -------
    np.ndarray
        LOS scalar field interpolated onto the full raster grid,
        shape ``(ny, nx)``.

    Raises
    ------
    ValueError
        If no valid GNSS LOS samples remain after masking.

    Examples
    --------
    >>> import geopandas as gpd
    >>> import numpy as np
    >>> from shapely.geometry import Point
    >>> from venti.gnss.los import project_to_los
    >>> stations = gpd.GeoDataFrame(
    ...     {"deast": [-14.0, -13.0, -14.5], "dnorth": [-3.0, -3.0, -2.5],
    ...      "dup": [-1.0, 0.5, -2.0]},
    ...     geometry=[Point(100, 100), Point(700, 400), Point(300, 500)],
    ... )
    >>> x, y = np.arange(0.0, 800, 10), np.arange(600.0, 0, -10)
    >>> los = [np.full((y.size, x.size), c) for c in (-0.6, -0.1, 0.78)]
    >>> field = project_to_los(*los, (x, y), stations, method="griddata")
    >>> field.shape
    (60, 80)

    """
    dtype = np.float32
    coords = grid_coordinates(grid)

    geo_x = np.asarray([p.x for p in gnss_gdf.geometry], dtype=dtype)
    geo_y = np.asarray([p.y for p in gnss_gdf.geometry], dtype=dtype)
    E = np.asarray(gnss_gdf["deast"].values, dtype=dtype)
    N = np.asarray(gnss_gdf["dnorth"].values, dtype=dtype)
    U = np.asarray(gnss_gdf["dup"].values, dtype=dtype)

    de_s, dn_s, du_s = _sample_los(
        coords, los_east, los_north, los_up, geo_x, geo_y, dtype
    )

    los_at_stations = (de_s * E + dn_s * N + du_s * U).astype(dtype)

    valid = np.isfinite(los_at_stations) & np.isfinite(geo_x) & np.isfinite(geo_y)
    gx = geo_x[valid]
    gy = geo_y[valid]
    zv = los_at_stations[valid]

    if zv.size == 0:
        msg = (
            "No valid GNSS LOS samples after masking — check CRS, coverage, and masks."
        )
        raise ValueError(msg)

    logger.info(
        "Interpolating LOS from %d GNSS stations using method='%s'", zv.size, method
    )

    if method == "rbf":
        return interpolate_rbf(coords, gx, gy, zv, function=rbf_function)
    if method == "griddata":
        return interpolate_griddata(coords, gx, gy, zv)

    msg = f"Unknown interpolation method '{method}'. Choose 'rbf' or 'griddata'."
    raise ValueError(msg)


def project_uncertainty_to_los(
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    grid: GridLike,
    gnss_gdf: gpd.GeoDataFrame,
    method: InterpolationMethod = "rbf",
    rbf_function: RbfFunction = "cubic",
) -> np.ndarray:
    """Project GNSS east/north/up uncertainties into an LOS uncertainty field.

    Companion to `project_to_los`: samples the LOS unit vectors at each GNSS
    station location, combines the east/north/up uncertainties in
    quadrature (assuming independent component errors) to get each
    station's own LOS rate uncertainty, then interpolates that scalar field
    onto the full raster grid the same way the LOS value itself is
    interpolated. Intended for `fit_windowed_surface`'s `gnss_los_std`:
    stations with a poorly-constrained rate (short record, few
    observations) should influence the windowed fit less than
    well-constrained ones.

    Parameters
    ----------
    los_east : np.ndarray
        2-D array of east LOS unit-vector components, shape ``(ny, nx)``.
    los_north : np.ndarray
        2-D array of north LOS unit-vector components.
    los_up : np.ndarray
        2-D array of up LOS unit-vector components.
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        Grid of the LOS arrays and the output: NetCDF, GeoTIFF, Dataset or
        ``(x, y)`` coordinate arrays, in the stations' CRS (see
        `venti.spatial.interpolation.grid_coordinates`).
    gnss_gdf : gpd.GeoDataFrame
        GNSS station data with a ``geometry`` column (Point, UTM) and
        columns ``'dsigma_e'``, ``'dsigma_n'``, ``'dsigma_u'`` containing
        the uncertainty of each ENU component.
    method : {'rbf', 'griddata'}, optional
        Spatial interpolation method, by default ``'rbf'``.
    rbf_function : str, optional
        RBF kernel passed to :func:`~venti.spatial.interpolation.interpolate_rbf`
        when ``method='rbf'``, by default ``'cubic'``.

    Returns
    -------
    np.ndarray
        LOS uncertainty field interpolated onto the full raster grid,
        shape ``(ny, nx)``, same units as `gnss_gdf`'s sigma columns.

    Raises
    ------
    ValueError
        If no valid GNSS LOS uncertainty samples remain after masking.

    """
    dtype = np.float32
    coords = grid_coordinates(grid)

    geo_x = np.asarray([p.x for p in gnss_gdf.geometry], dtype=dtype)
    geo_y = np.asarray([p.y for p in gnss_gdf.geometry], dtype=dtype)
    SE = np.asarray(gnss_gdf["dsigma_e"].values, dtype=dtype)
    SN = np.asarray(gnss_gdf["dsigma_n"].values, dtype=dtype)
    SU = np.asarray(gnss_gdf["dsigma_u"].values, dtype=dtype)

    de_s, dn_s, du_s = _sample_los(
        coords, los_east, los_north, los_up, geo_x, geo_y, dtype
    )

    # Independent-error propagation through the same LOS dot product
    # project_to_los uses for the value itself.
    los_std_at_stations = np.sqrt(
        (de_s * SE) ** 2 + (dn_s * SN) ** 2 + (du_s * SU) ** 2
    ).astype(dtype)

    valid = np.isfinite(los_std_at_stations) & np.isfinite(geo_x) & np.isfinite(geo_y)
    gx = geo_x[valid]
    gy = geo_y[valid]
    zv = los_std_at_stations[valid]

    if zv.size == 0:
        msg = (
            "No valid GNSS LOS uncertainty samples after masking — "
            "check CRS, coverage, and masks."
        )
        raise ValueError(msg)

    logger.info(
        "Interpolating LOS uncertainty from %d GNSS stations using method='%s'",
        zv.size,
        method,
    )

    if method == "rbf":
        return interpolate_rbf(coords, gx, gy, zv, function=rbf_function)
    if method == "griddata":
        return interpolate_griddata(coords, gx, gy, zv)

    msg = f"Unknown interpolation method '{method}'. Choose 'rbf' or 'griddata'."
    raise ValueError(msg)
