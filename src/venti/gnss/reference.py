# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""GNSSReference: high-level interface for GNSS-based InSAR calibration data."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    import geopandas as gpd

    from ..spatial.interpolation import GridLike

from .los import (
    InterpolationMethod,
    RbfFunction,
    project_to_los,
    project_uncertainty_to_los,
)
from .unr import (
    calculate_station_velocity,
    download_grid_lookup,
    download_station,
    find_stations_in_bounds,
    read_epoch_displacements,
    read_station_rate,
)

logger = logging.getLogger(__name__)


@dataclass
class GNSSReference:
    """GNSS grid reference data for a geographic region.

    Downloads UNR GNSS grid timeseries for stations within the specified
    bounds and provides methods for projecting GNSS observations into the
    InSAR line-of-sight (LOS) direction.

    Parameters
    ----------
    bounds : tuple of float
        Bounding box as ``(south, north, west, east)`` in the UTM CRS.
    output_dir : Path
        Directory for downloaded GNSS files.
    reference_frame : str, optional
        GNSS reference frame, ``'IGS20'`` or ``'IGS14'``, by default ``'IGS20'``.
    utm_epsg : int or None, optional
        EPSG code of the UTM projection. Required before calling
        :meth:`download_stations`.
    grid_type : str, optional
        Which UNR grid product :meth:`download_stations` fetches:
        ``'constant'`` (precomputed linear rates, IGS20 only) or
        ``'variable'`` (per-epoch positions), by default ``'constant'``.
        Determines which downstream computation is valid:
        :meth:`compute_velocity_los` for ``'constant'``,
        :meth:`compute_displacement_los` for ``'variable'``.

    Attributes
    ----------
    station_files : list of Path
        Paths to downloaded station files (set after :meth:`download_stations`).
    station_gdf : gpd.GeoDataFrame or None
        Filtered station GeoDataFrame in UTM CRS (set after :meth:`download_stations`).

    Examples
    --------
    Download stations and compute a constant LOS velocity field::

        gnss = GNSSReference(
            bounds=(3800000, 3900000, 400000, 500000),
            output_dir=Path('output/GNSS'),
            utm_epsg=32611,
        )
        n = gnss.download_stations()
        gnss_los = gnss.compute_velocity_los(
            los_east, los_north, los_up, 'displacement.nc'
        )

    """

    bounds: tuple[float, float, float, float]
    output_dir: Path
    reference_frame: str = "IGS20"
    utm_epsg: int | None = None
    grid_type: str = "constant"

    station_files: list[Path] = field(default_factory=list, init=False)
    station_gdf: gpd.GeoDataFrame | None = field(default=None, init=False)

    def download_stations(self) -> int:
        """Download all GNSS stations within the configured bounds.

        Returns
        -------
        int
            Number of stations downloaded (or already cached).

        Raises
        ------
        ValueError
            If :attr:`utm_epsg` is not set.

        """
        if self.utm_epsg is None:
            msg = "utm_epsg must be set before calling download_stations"
            raise ValueError(msg)

        self.output_dir.mkdir(parents=True, exist_ok=True)

        lookup_path = self.output_dir / "grid_latlon_lookup.txt"
        if not lookup_path.exists():
            lookup_path = download_grid_lookup(self.output_dir, self.reference_frame)

        self.station_gdf = find_stations_in_bounds(
            lookup_path, self.bounds, self.utm_epsg
        )
        logger.info("Found %d GNSS stations in bounds", len(self.station_gdf))

        self.station_files = []
        for station_id in self.station_gdf.index:
            path = download_station(
                station_id,
                self.output_dir,
                self.reference_frame,
                grid_type=self.grid_type,
            )
            self.station_files.append(path)

        logger.info("Downloaded %d station files", len(self.station_files))
        return len(self.station_files)

    def _build_velocity_gdf(self, start_year: float = 2014.0):
        """Compute velocity GeoDataFrame for all stations.

        For ``grid_type='constant'``, rates and their uncertainties are read
        directly from the precomputed constant-grid files (`start_year` is
        not applicable and is ignored — see :func:`~.unr.read_station_rate`).
        For ``grid_type='variable'``, rates are estimated by weighted linear
        regression over the per-epoch positions (see
        :func:`~.unr.calculate_station_velocity`).

        Parameters
        ----------
        start_year : float, optional
            Exclude observations before this decimal year, by default
            ``2014.0``. Only applies when `grid_type` is ``'variable'``.

        Returns
        -------
        gpd.GeoDataFrame
            Velocity GeoDataFrame with columns ``deast``, ``dnorth``, ``dup``,
            ``dsigma_e``, ``dsigma_n``, ``dsigma_u``, and ``geometry``.

        """
        import geopandas as gpd
        import pandas as pd

        if self.grid_type == "constant":
            logger.info(
                "grid_type='constant': reading precomputed rates directly; "
                "start_year=%.4f is not applied",
                start_year,
            )

        rows = []
        for path in self.station_files:
            station_id = int(path.name.split("_")[0])
            if self.grid_type == "constant":
                ve, vn, vu, se, sn, su = read_station_rate(path)
            else:
                ve, vn, vu, se, sn, su = calculate_station_velocity(path, start_year)
            rows.append(
                {
                    "id": station_id,
                    "deast": ve,
                    "dnorth": vn,
                    "dup": vu,
                    "dsigma_e": se,
                    "dsigma_n": sn,
                    "dsigma_u": su,
                }
            )

        assert self.station_gdf is not None, "Call download_stations() first"
        df = pd.DataFrame(rows).set_index("id")
        return gpd.GeoDataFrame(df.join(self.station_gdf[["geometry"]], how="inner"))

    def compute_velocity_los(
        self,
        los_east: np.ndarray,
        los_north: np.ndarray,
        los_up: np.ndarray,
        grid: GridLike,
        start_year: float = 2014.0,
        method: InterpolationMethod = "rbf",
        rbf_function: RbfFunction = "cubic",
    ) -> np.ndarray:
        """Compute GNSS velocity projected into the InSAR LOS direction.

        Parameters
        ----------
        los_east : np.ndarray
            2-D east LOS unit-vector component, shape ``(ny, nx)``.
        los_north : np.ndarray
            2-D north LOS unit-vector component.
        los_up : np.ndarray
            2-D up LOS unit-vector component.
        grid : str, Path, xr.Dataset or tuple of np.ndarray
            Grid of the LOS arrays and output (see `compute_gnss_los`).
        start_year : float, optional
            Earliest observation year used in velocity estimation,
            by default ``2014.0``.
        method : {'rbf', 'griddata'}, optional
            Spatial interpolation method, by default ``'rbf'``.
        rbf_function : str, optional
            RBF basis function when ``method='rbf'``, by default ``'cubic'``.

        Returns
        -------
        np.ndarray
            GNSS LOS velocity field, shape ``(ny, nx)``, in mm/yr.

        """
        assert self.station_files, "Call download_stations() first"
        velocity_gdf = self._build_velocity_gdf(start_year)
        return project_to_los(
            los_east,
            los_north,
            los_up,
            grid,
            velocity_gdf,
            method=method,
            rbf_function=rbf_function,
        )

    def compute_velocity_los_std(
        self,
        los_east: np.ndarray,
        los_north: np.ndarray,
        los_up: np.ndarray,
        grid: GridLike,
        start_year: float = 2014.0,
        method: InterpolationMethod = "rbf",
        rbf_function: RbfFunction = "cubic",
    ) -> np.ndarray:
        """Compute GNSS velocity uncertainty projected into the InSAR LOS direction.

        Companion to `compute_velocity_los`, for `fit_windowed_surface`'s
        `gnss_los_std`. Uses the same station set and the same interpolation
        method, so the value and uncertainty fields are built consistently
        with each other.

        Parameters
        ----------
        los_east : np.ndarray
            2-D east LOS unit-vector component, shape ``(ny, nx)``.
        los_north : np.ndarray
            2-D north LOS unit-vector component.
        los_up : np.ndarray
            2-D up LOS unit-vector component.
        grid : str, Path, xr.Dataset or tuple of np.ndarray
            Grid of the LOS arrays and output (see `compute_gnss_los`).
        start_year : float, optional
            Earliest observation year used in velocity estimation,
            by default ``2014.0``.
        method : {'rbf', 'griddata'}, optional
            Spatial interpolation method, by default ``'rbf'``.
        rbf_function : str, optional
            RBF basis function when ``method='rbf'``, by default ``'cubic'``.

        Returns
        -------
        np.ndarray
            GNSS LOS velocity uncertainty field, shape ``(ny, nx)``, in mm/yr.

        """
        assert self.station_files, "Call download_stations() first"
        velocity_gdf = self._build_velocity_gdf(start_year)
        return project_uncertainty_to_los(
            los_east,
            los_north,
            los_up,
            grid,
            velocity_gdf,
            method=method,
            rbf_function=rbf_function,
        )

    def compute_displacement_los(
        self,
        ref_date: float,
        sec_date: float,
        los_east: np.ndarray,
        los_north: np.ndarray,
        los_up: np.ndarray,
        grid: GridLike,
        method: InterpolationMethod = "rbf",
        rbf_function: RbfFunction = "cubic",
    ) -> np.ndarray:
        """Compute epoch-specific GNSS displacement projected into LOS.

        Parameters
        ----------
        ref_date : float
            Reference epoch as decimal year.
        sec_date : float
            Secondary epoch as decimal year.
        los_east : np.ndarray
            2-D east LOS unit-vector component.
        los_north : np.ndarray
            2-D north LOS unit-vector component.
        los_up : np.ndarray
            2-D up LOS unit-vector component.
        grid : str, Path, xr.Dataset or tuple of np.ndarray
            Grid of the LOS arrays and output (see `compute_gnss_los`).
        method : {'rbf', 'griddata'}, optional
            Spatial interpolation method, by default ``'rbf'``.
        rbf_function : str, optional
            RBF basis function when ``method='rbf'``, by default ``'cubic'``.

        Returns
        -------
        np.ndarray
            GNSS LOS displacement field, shape ``(ny, nx)``.

        """
        assert self.station_files, "Call download_stations() first"
        assert self.station_gdf is not None, "Call download_stations() first"

        disp_gdf = read_epoch_displacements(
            self.station_files, ref_date, sec_date, self.station_gdf
        )
        return project_to_los(
            los_east,
            los_north,
            los_up,
            grid,
            disp_gdf,
            method=method,
            rbf_function=rbf_function,
        )

    def compute_displacement_los_std(
        self,
        ref_date: float,
        sec_date: float,
        los_east: np.ndarray,
        los_north: np.ndarray,
        los_up: np.ndarray,
        grid: GridLike,
        method: InterpolationMethod = "rbf",
        rbf_function: RbfFunction = "cubic",
    ) -> np.ndarray:
        """Compute epoch-specific GNSS displacement uncertainty projected into LOS.

        Companion to `compute_displacement_los`, built from the same stations
        and epochs: each station's ``ref`` and ``sec`` position sigmas are
        combined in quadrature per component, then projected like
        `compute_velocity_los_std`.

        Parameters
        ----------
        ref_date : float
            Reference epoch as decimal year.
        sec_date : float
            Secondary epoch as decimal year.
        los_east : np.ndarray
            2-D east LOS unit-vector component.
        los_north : np.ndarray
            2-D north LOS unit-vector component.
        los_up : np.ndarray
            2-D up LOS unit-vector component.
        grid : str, Path, xr.Dataset or tuple of np.ndarray
            Grid of the LOS arrays and output (see `compute_gnss_los`).
        method : {'rbf', 'griddata'}, optional
            Spatial interpolation method, by default ``'rbf'``.
        rbf_function : str, optional
            RBF basis function when ``method='rbf'``, by default ``'cubic'``.

        Returns
        -------
        np.ndarray
            GNSS LOS displacement uncertainty field, shape ``(ny, nx)``, in mm.

        """
        assert self.station_files, "Call download_stations() first"
        assert self.station_gdf is not None, "Call download_stations() first"

        disp_gdf = read_epoch_displacements(
            self.station_files, ref_date, sec_date, self.station_gdf
        )
        return project_uncertainty_to_los(
            los_east,
            los_north,
            los_up,
            grid,
            disp_gdf,
            method=method,
            rbf_function=rbf_function,
        )


def _load_cache(cache_file: Path, expected_shape: tuple[int, ...]) -> np.ndarray:
    """Load a ``.npy`` GNSS cache, refusing one built for a different grid.

    Cache files are keyed only by directory and date, so one written for a
    different LOS raster would otherwise be silently reused.

    Raises
    ------
    ValueError
        If the cached array's shape differs from `expected_shape`.

    """
    cached = np.load(cache_file)
    if cached.shape != expected_shape:
        msg = (
            f"Stale GNSS cache {cache_file}: shape {cached.shape} does not match "
            f"the LOS grid {expected_shape}. Delete it or set recompute_gnss=True."
        )
        raise ValueError(msg)
    return cached


def _cached(
    cache_dir: Path | None,
    name: str,
    shape: tuple[int, ...],
    recompute: bool,
    compute: Callable[[], np.ndarray],
    what: str,
) -> np.ndarray:
    """Load ``cache_dir/name`` if present, else compute (and save if cached)."""
    if cache_dir is None:
        return compute()
    cache_file = Path(cache_dir) / name
    if cache_file.exists() and not recompute:
        logger.info("Loading cached %s from %s", what, cache_file)
        return _load_cache(cache_file, shape)
    result = compute()
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache_file, result)
    logger.info("Saved %s to %s", what, cache_file)
    return result


def _require_dates(ref_date: float | None, sec_date: float | None) -> None:
    if ref_date is None or sec_date is None:
        msg = "ref_date and sec_date are required for grid_type='variable'"
        raise ValueError(msg)


def compute_gnss_los(
    gnss_ref: GNSSReference,
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    grid: GridLike,
    cache_dir: Path | None = None,
    ref_date: float | None = None,
    sec_date: float | None = None,
    recompute: bool = False,
) -> np.ndarray:
    """Compute the GNSS LOS displacement (or velocity) field on a grid.

    The mode follows the product `gnss_ref` downloaded (`gnss_ref.grid_type`):

    - ``'constant'``: the LOS velocity field (mm/yr) from UNR's rate grid,
      times ``sec_date - ref_date`` (mm); the velocity itself if no dates are
      given. The velocity depends only on geometry and is cached once as
      ``gnss_los_velocity.npy``.
    - ``'variable'``: the LOS displacement between the positions nearest
      the two dates (mm), cached per pair as
      ``gnss_los_disp_<ref>_<sec>.npy``.

    Both are ``sec - ref``, in the grid's own frame (IGS20), plate motion
    included; separating vertical from horizontal motion is left to the LOS
    decomposition.

    Parameters
    ----------
    gnss_ref : GNSSReference
        GNSS reference with stations already downloaded.
    los_east : np.ndarray
        2-D east LOS unit-vector component, shape ``(ny, nx)``.
    los_north : np.ndarray
        2-D north LOS unit-vector component.
    los_up : np.ndarray
        2-D up LOS unit-vector component.
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        Grid of the LOS arrays and output (NetCDF, GeoTIFF, Dataset or
        ``(x, y)`` arrays; see `venti.spatial.interpolation.grid_coordinates`).
    cache_dir : Path, optional
        Directory for ``.npy`` caches. ``None`` (default) disables caching.
    ref_date : float, optional
        Reference epoch as decimal year. Required for ``'variable'``.
    sec_date : float, optional
        Secondary epoch as decimal year. Required for ``'variable'``.
    recompute : bool, optional
        Ignore existing cache files, by default ``False``.

    Returns
    -------
    np.ndarray
        GNSS LOS field in mm (mm/yr for ``'constant'`` without dates).

    Raises
    ------
    ValueError
        If dates are missing for ``'variable'``, or a cache does not match
        the grid shape.

    Examples
    --------
    ::

        gnss = GNSSReference(bounds, Path("gnss"), utm_epsg=32618)
        gnss.download_stations()
        gnss_los_mm = compute_gnss_los(
            gnss, los_e, los_n, los_u, "product.nc",
            ref_date=2020.0, sec_date=2020.5,
        )

    """
    if gnss_ref.grid_type == "constant":
        velocity = _cached(
            cache_dir,
            "gnss_los_velocity.npy",
            los_east.shape,
            recompute,
            lambda: gnss_ref.compute_velocity_los(
                los_east=los_east,
                los_north=los_north,
                los_up=los_up,
                grid=grid,
                # Govorcin et al. (2025) use linear sampling to avoid RBF
                # overshoot; against independent stations the two match in
                # bias/RMSE, and RBF is smoother between lattice points.
                method="rbf",
            ),
            "GNSS LOS velocity",
        )
        if ref_date is None or sec_date is None:
            return velocity
        # sec_date > ref_date (OPERA convention), so the interval is positive.
        return velocity * (sec_date - ref_date)

    _require_dates(ref_date, sec_date)
    assert ref_date is not None
    assert sec_date is not None
    return _cached(
        cache_dir,
        f"gnss_los_disp_{ref_date:.4f}_{sec_date:.4f}.npy",
        los_east.shape,
        recompute,
        lambda: gnss_ref.compute_displacement_los(
            ref_date=ref_date,
            sec_date=sec_date,
            los_east=los_east,
            los_north=los_north,
            los_up=los_up,
            grid=grid,
            method="rbf",
        ),
        "GNSS LOS displacement",
    )


def compute_gnss_los_std(
    gnss_ref: GNSSReference,
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    grid: GridLike,
    cache_dir: Path | None = None,
    ref_date: float | None = None,
    sec_date: float | None = None,
    recompute: bool = False,
) -> np.ndarray:
    """Compute the GNSS LOS uncertainty field, for weighted fitting.

    Companion to `compute_gnss_los` with the same modes, grid and caching:
    the rate uncertainty times the interval for ``'constant'`` (cache
    ``gnss_los_velocity_std.npy``), the combined ``ref``/``sec`` position
    uncertainty for ``'variable'`` (cache ``gnss_los_disp_std_<ref>_<sec>.npy``).

    Parameters
    ----------
    gnss_ref : GNSSReference
        GNSS reference with stations already downloaded.
    los_east : np.ndarray
        2-D east LOS unit-vector component, shape ``(ny, nx)``.
    los_north : np.ndarray
        2-D north LOS unit-vector component.
    los_up : np.ndarray
        2-D up LOS unit-vector component.
    grid : str, Path, xr.Dataset or tuple of np.ndarray
        Grid of the LOS arrays and output (see `compute_gnss_los`).
    cache_dir : Path, optional
        Directory for ``.npy`` caches. ``None`` (default) disables caching.
    ref_date : float, optional
        Reference epoch as decimal year. Required for ``'variable'``.
    sec_date : float, optional
        Secondary epoch as decimal year. Required for ``'variable'``.
    recompute : bool, optional
        Ignore existing cache files, by default ``False``.

    Returns
    -------
    np.ndarray
        GNSS LOS uncertainty in mm (mm/yr for ``'constant'`` without dates).

    Raises
    ------
    ValueError
        If dates are missing for ``'variable'``, or a cache does not match
        the grid shape.

    """
    if gnss_ref.grid_type == "constant":
        velocity_std = _cached(
            cache_dir,
            "gnss_los_velocity_std.npy",
            los_east.shape,
            recompute,
            lambda: gnss_ref.compute_velocity_los_std(
                los_east=los_east,
                los_north=los_north,
                los_up=los_up,
                grid=grid,
                method="rbf",
            ),
            "GNSS LOS velocity uncertainty",
        )
        if ref_date is None or sec_date is None:
            return velocity_std
        return velocity_std * (sec_date - ref_date)

    _require_dates(ref_date, sec_date)
    assert ref_date is not None
    assert sec_date is not None
    return _cached(
        cache_dir,
        f"gnss_los_disp_std_{ref_date:.4f}_{sec_date:.4f}.npy",
        los_east.shape,
        recompute,
        lambda: gnss_ref.compute_displacement_los_std(
            ref_date=ref_date,
            sec_date=sec_date,
            los_east=los_east,
            los_north=los_north,
            los_up=los_up,
            grid=grid,
            method="rbf",
        ),
        "GNSS LOS uncertainty",
    )
