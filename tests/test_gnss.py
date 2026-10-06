# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Unit tests for the gnss subpackage.

Covers pure-computation functions only — no network calls.
Download helpers (download_grid_lookup, download_station) are excluded
because they require internet access.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

try:
    import geopandas as gpd
    from shapely.geometry import Point

    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False

try:
    import xarray as xr

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False

from venti.gnss.los import project_to_los
from venti.gnss.unr import (
    calculate_station_velocity,
    download_station,
    find_stations_in_bounds,
    read_epoch_displacements,
    read_station_rate,
)
from venti.spatial.interpolation import (
    _regular_grid_interpolator,
    _sample_on_points,
    interpolate_griddata,
    interpolate_rbf,
)


# Helpers shared across tests
def _write_tenv8(path: Path, years: np.ndarray, east, north, up, sigma=1e-3) -> None:
    """Write a minimal UNR .tenv8-style whitespace-delimited station file."""
    rows = np.column_stack(
        [
            years,
            east,
            north,
            up,
            np.full_like(years, sigma),
            np.full_like(years, sigma),
            np.full_like(years, sigma),
        ]
    )
    np.savetxt(path, rows, fmt="%.6f")


def _make_netcdf(
    path: Path, nx: int = 20, ny: int = 15
) -> tuple[np.ndarray, np.ndarray]:
    """Create a minimal NetCDF with x/y UTM coords; return (x, y) arrays."""
    x = np.linspace(400_000, 500_000, nx, dtype=np.float32)
    y = np.linspace(3_800_000, 3_900_000, ny, dtype=np.float32)
    data = np.ones((ny, nx), dtype=np.float32)
    ds = xr.Dataset({"displacement": (["y", "x"], data)}, coords={"x": x, "y": y})
    ds.to_netcdf(path)
    return x, y


def _make_gnss_gdf(
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    deast: np.ndarray,
    dnorth: np.ndarray,
    dup: np.ndarray,
) -> gpd.GeoDataFrame:
    """Build a minimal GNSS GeoDataFrame in the same CRS as _make_netcdf."""
    return gpd.GeoDataFrame(
        {
            "deast": deast,
            "dnorth": dnorth,
            "dup": dup,
        },
        geometry=[Point(x, y) for x, y in zip(x_coords, y_coords, strict=False)],
        crs="EPSG:32611",
    )


# venti.gnss.unr
@pytest.mark.skipif(not HAS_GEOPANDAS, reason="geopandas not installed")
class TestFindStationsInBounds:
    """Tests for find_stations_in_bounds."""

    def _write_lookup(self, path: Path, records: list[tuple]) -> None:
        """Write a minimal grid_latlon_lookup.txt (id lon lat)."""
        lines = [f"{id_} {lon} {lat}" for id_, lon, lat in records]
        path.write_text("\n".join(lines))

    def test_stations_within_bounds_returned(self):
        with tempfile.TemporaryDirectory() as tmp:
            lookup = Path(tmp) / "grid_latlon_lookup.txt"
            # UTM zone 11N: approx lon -120..-114, lat 34..36
            self._write_lookup(
                lookup,
                [
                    (1, -118.0, 35.0),  # inside: easting~409k, northing~3870k
                    (2, -117.5, 34.5),  # inside: easting~455k, northing~3814k
                    (3, -110.0, 35.0),  # outside (east)
                ],
            )
            # bounds in UTM zone 11N metres (S, N, W, E)
            bounds = (3_800_000, 3_900_000, 350_000, 500_000)
            result = find_stations_in_bounds(lookup, bounds, utm_epsg=32611)
            assert len(result) == 2
            assert 3 not in result.index

    def test_empty_result_when_no_stations_in_bounds(self):
        with tempfile.TemporaryDirectory() as tmp:
            lookup = Path(tmp) / "grid_latlon_lookup.txt"
            self._write_lookup(lookup, [(1, -90.0, 35.0)])
            bounds = (3_800_000, 3_900_000, 350_000, 500_000)
            result = find_stations_in_bounds(lookup, bounds, utm_epsg=32611)
            assert len(result) == 0

    def test_padding_extends_search_area(self):
        with tempfile.TemporaryDirectory() as tmp:
            lookup = Path(tmp) / "grid_latlon_lookup.txt"
            self._write_lookup(lookup, [(1, -118.0, 35.0)])
            bounds = (3_800_000, 3_900_000, 350_000, 500_000)
            no_pad = find_stations_in_bounds(lookup, bounds, utm_epsg=32611)
            with_pad = find_stations_in_bounds(
                lookup, bounds, utm_epsg=32611, padding=200_000
            )
            assert len(with_pad) >= len(no_pad)

    def test_output_indexed_by_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            lookup = Path(tmp) / "grid_latlon_lookup.txt"
            self._write_lookup(lookup, [(42, -118.0, 35.0)])
            bounds = (3_800_000, 3_900_000, 350_000, 500_000)
            result = find_stations_in_bounds(lookup, bounds, utm_epsg=32611)
            assert 42 in result.index


class TestCalculateStationVelocity:
    """Tests for calculate_station_velocity — numpy only, no optional deps."""

    def test_known_linear_trend(self):
        """Velocity recovered from a perfect linear trend should match the slope."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "000001_IGS20.tenv8"
            years = np.linspace(2014.0, 2024.0, 50)
            ve_true, vn_true, vu_true = 0.010, -0.005, 0.002  # m/yr
            _write_tenv8(
                path,
                years,
                ve_true * years + 0.1,
                vn_true * years + 0.05,
                vu_true * years - 0.2,
                sigma=1e-4,
            )
            ve, vn, vu, *_ = calculate_station_velocity(path, start_year=2014.0)
            assert abs(ve - ve_true) < 1e-4
            assert abs(vn - vn_true) < 1e-4
            assert abs(vu - vu_true) < 1e-4

    def test_start_year_filter(self):
        """Observations before start_year must be excluded."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "000002_IGS20.tenv8"
            years_early = np.linspace(2010.0, 2013.9, 20)
            years_late = np.linspace(2014.0, 2024.0, 50)
            years = np.concatenate([years_early, years_late])
            east = np.concatenate([np.zeros(20), 0.01 * years_late])
            _write_tenv8(path, years, east, np.zeros_like(years), np.zeros_like(years))
            ve, *_ = calculate_station_velocity(path, start_year=2014.0)
            assert ve > 0.005

    def test_returns_six_floats(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "000003_IGS20.tenv8"
            years = np.linspace(2014.0, 2024.0, 30)
            _write_tenv8(path, years, years * 0.01, years * 0.005, years * 0.002)
            result = calculate_station_velocity(path)
            assert len(result) == 6
            assert all(isinstance(v, float) for v in result)


@pytest.mark.skipif(not HAS_GEOPANDAS, reason="geopandas not installed")
class TestReadEpochDisplacements:
    """Tests for read_epoch_displacements."""

    def _make_station_gdf(self, station_ids: list[int]) -> gpd.GeoDataFrame:
        return gpd.GeoDataFrame(
            index=station_ids,
            geometry=[Point(400_000 + i * 10_000, 3_850_000) for i in station_ids],
            crs="EPSG:32611",
        )

    def test_displacement_is_sec_minus_ref(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "000001_IGS20.tenv8"
            years = np.array([2020.0, 2020.5, 2021.0])
            east = np.array([0.00, 0.01, 0.02])
            _write_tenv8(path, years, east, east * 0.5, east * 0.1)

            gdf = self._make_station_gdf([1])
            result = read_epoch_displacements([path], 2020.0, 2021.0, gdf)

            # Same sign as OPERA displacement: position at sec minus at ref.
            assert abs(result.loc[1, "deast"] - (0.02 - 0.00)) < 1e-6

    def test_matches_constant_rate_times_interval(self):
        # The same linear motion must give the same GNSS displacement whether
        # read as per-epoch positions ('variable') or as a rate ('constant').
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "000001_IGS20.tenv8"
            years = 2016.0 + np.arange(0, 6 * 365) / 365.25
            ve, vn, vu = -14.0, -3.0, -1.5
            _write_tenv8(
                path,
                years,
                ve * (years - 2016),
                vn * (years - 2016),
                vu * (years - 2016),
            )
            ref_date, sec_date = years[100], years[900]

            gdf = self._make_station_gdf([1])
            disp = read_epoch_displacements([path], ref_date, sec_date, gdf).loc[1]
            rate = read_station_rate(path)

            for got, v in zip(
                (disp["deast"], disp["dnorth"], disp["dup"]), rate[:3], strict=False
            ):
                assert got == pytest.approx(v * (sec_date - ref_date), abs=1e-6)

    def test_station_without_nearby_epoch_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            covered = Path(tmp) / "000001_IGS20.tenv8"
            late = Path(tmp) / "000002_IGS20.tenv8"
            _write_tenv8(
                covered,
                np.array([2020.0, 2021.0]),
                np.zeros(2),
                np.zeros(2),
                np.zeros(2),
            )
            _write_tenv8(
                late, np.array([2020.8, 2021.0]), np.zeros(2), np.zeros(2), np.zeros(2)
            )

            gdf = self._make_station_gdf([1, 2])
            result = read_epoch_displacements([covered, late], 2020.0, 2021.0, gdf)
            assert list(result.index) == [1]

            with pytest.raises(ValueError, match="No GNSS station covers"):
                read_epoch_displacements([late], 2020.0, 2021.0, gdf)

    def test_output_has_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "000001_IGS20.tenv8"
            years = np.array([2020.0, 2021.0])
            _write_tenv8(path, years, years * 0.01, years * 0.01, years * 0.01)
            gdf = self._make_station_gdf([1])
            result = read_epoch_displacements([path], 2020.0, 2021.0, gdf)
            assert "geometry" in result.columns

    def test_multiple_stations(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for i in [1, 2]:
                p = Path(tmp) / f"{i:06d}_IGS20.tenv8"
                years = np.array([2020.0, 2021.0])
                _write_tenv8(p, years, years * 0.01 * i, years * 0.01, years * 0.01)
                paths.append(p)

            gdf = self._make_station_gdf([1, 2])
            result = read_epoch_displacements(paths, 2020.0, 2021.0, gdf)
            assert len(result) == 2


def test_station_cache_is_separate_per_grid_type(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from venti.gnss import unr

    requested = []

    def fake_get(url, timeout):
        requested.append(url)
        kind = "variable" if "time_variable" in url else "constant"
        return SimpleNamespace(ok=True, status_code=200, text=kind)

    monkeypatch.setattr(unr.requests, "get", fake_get)
    constant = download_station(1, tmp_path, "IGS20", grid_type="constant")
    variable = download_station(1, tmp_path, "IGS20", grid_type="variable")

    assert constant != variable
    assert constant.read_text() == "constant"
    assert variable.read_text() == "variable"
    # Both products share the remote filename.
    assert all(url.endswith("/000001_IGS20.tenv8") for url in requested)


def test_variable_grid_gnss_uncertainty_is_computed_per_epoch(tmp_path):
    from types import SimpleNamespace

    from venti.gnss.reference import compute_gnss_los_std

    calls = []

    def fake_std(**kwargs):
        calls.append((kwargs["ref_date"], kwargs["sec_date"]))
        return np.full((3, 4), 2.0)

    gnss_ref = SimpleNamespace(
        grid_type="variable", compute_displacement_los_std=fake_std
    )
    kwargs = {
        "gnss_ref": gnss_ref,
        "los_east": np.zeros((3, 4)),
        "los_north": np.zeros((3, 4)),
        "los_up": np.zeros((3, 4)),
        "grid": tmp_path / "disp.nc",
        "cache_dir": tmp_path,
        "ref_date": 2020.0,
        "sec_date": 2020.5,
    }
    first = compute_gnss_los_std(**kwargs)
    cached = compute_gnss_los_std(**kwargs)

    np.testing.assert_array_equal(first, np.full((3, 4), 2.0))
    np.testing.assert_array_equal(cached, first)
    assert calls == [(2020.0, 2020.5)]


# venti.gnss.los
@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
class TestRegularGridInterpolatorFromDs:
    """Tests for _regular_grid_interpolator."""

    def test_returns_callable(self):
        x = np.linspace(0, 10, 5)
        y = np.linspace(0, 10, 4)
        arr = np.ones((4, 5), dtype=np.float32)
        ds = xr.Dataset({"v": (["y", "x"], arr)}, coords={"x": x, "y": y})
        interp = _regular_grid_interpolator(ds, arr)
        pts = np.column_stack([[5.0], [5.0]])
        val = interp(pts)
        assert np.isfinite(val[0])

    def test_handles_decreasing_y(self):
        """y axis decreasing (north-up raster) must be flipped internally."""
        x = np.linspace(0, 10, 5)
        y = np.linspace(10, 0, 4)  # decreasing
        arr = np.ones((4, 5), dtype=np.float32)
        ds = xr.Dataset({"v": (["y", "x"], arr)}, coords={"x": x, "y": y})
        interp = _regular_grid_interpolator(ds, arr)
        pts = np.column_stack([[5.0], [5.0]])
        val = interp(pts)
        assert np.isfinite(val[0])


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
class TestSampleOnPoints:
    """Tests for _sample_on_points."""

    def test_sample_at_grid_node_returns_exact_value(self):
        x = np.array([0.0, 10.0, 20.0])
        y = np.array([0.0, 10.0])
        arr = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32)
        ds = xr.Dataset({"v": (["y", "x"], arr)}, coords={"x": x, "y": y})
        vals = _sample_on_points(ds, arr, pts_x=np.array([10.0]), pts_y=np.array([0.0]))
        assert abs(vals[0] - 2.0) < 1e-5

    def test_out_of_bounds_returns_nan(self):
        x = np.array([0.0, 10.0])
        y = np.array([0.0, 10.0])
        arr = np.ones((2, 2), dtype=np.float32)
        ds = xr.Dataset({"v": (["y", "x"], arr)}, coords={"x": x, "y": y})
        vals = _sample_on_points(
            ds, arr, pts_x=np.array([999.0]), pts_y=np.array([999.0])
        )
        assert np.isnan(vals[0])


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
class TestInterpolateRbf:
    """Tests for interpolate_rbf."""

    def test_output_shape_matches_netcdf_grid(self):
        with tempfile.TemporaryDirectory() as tmp:
            nc = Path(tmp) / "grid.nc"
            x, y = _make_netcdf(nc, nx=10, ny=8)
            gx = np.array([x[2], x[5], x[8]], dtype=np.float32)
            gy = np.array([y[1], y[4], y[6]], dtype=np.float32)
            zv = np.array([1.0, 2.0, 3.0], dtype=np.float32)
            result = interpolate_rbf(nc, gx, gy, zv, function="linear")
            assert result.shape == (8, 10)

    def test_constant_field_recovered(self):
        """RBF interpolation of a constant field should return that constant."""
        with tempfile.TemporaryDirectory() as tmp:
            nc = Path(tmp) / "grid.nc"
            x, y = _make_netcdf(nc, nx=10, ny=8)
            gx = np.array([x[0], x[-1], x[0], x[-1], x[5]], dtype=np.float32)
            gy = np.array([y[0], y[0], y[-1], y[-1], y[4]], dtype=np.float32)
            zv = np.full(5, 5.0, dtype=np.float32)
            result = interpolate_rbf(nc, gx, gy, zv, function="linear")
            interior = result[2:-2, 2:-2]
            assert np.allclose(interior, 5.0, atol=0.5)


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
class TestInterpolateGriddata:
    """Tests for interpolate_griddata."""

    def test_output_shape(self):
        with tempfile.TemporaryDirectory() as tmp:
            nc = Path(tmp) / "grid.nc"
            x, y = _make_netcdf(nc, nx=10, ny=8)
            gx = np.array([x[1], x[5], x[8]], dtype=np.float32)
            gy = np.array([y[1], y[4], y[6]], dtype=np.float32)
            zv = np.array([0.1, 0.2, 0.3], dtype=np.float32)
            result = interpolate_griddata(nc, gx, gy, zv, method="nearest")
            assert result.shape == (8, 10)


@pytest.mark.skipif(not HAS_GEOPANDAS, reason="geopandas not installed")
@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
class TestProjectToLos:
    """Tests for project_to_los."""

    def test_output_shape(self):
        with tempfile.TemporaryDirectory() as tmp:
            nc = Path(tmp) / "grid.nc"
            x, y = _make_netcdf(nc, nx=12, ny=10)

            ny, nx = 10, 12
            los_e = np.full((ny, nx), 0.6, dtype=np.float32)
            los_n = np.full((ny, nx), 0.1, dtype=np.float32)
            los_u = np.full((ny, nx), 0.8, dtype=np.float32)

            gdf = _make_gnss_gdf(
                x_coords=np.array([x[2], x[6], x[10]], dtype=np.float32),
                y_coords=np.array([y[3], y[6], y[8]], dtype=np.float32),
                deast=np.array([0.01, 0.02, 0.01]),
                dnorth=np.array([0.005, 0.005, 0.005]),
                dup=np.array([0.002, 0.002, 0.002]),
            )

            result = project_to_los(los_e, los_n, los_u, nc, gdf, method="griddata")
            assert result.shape == (ny, nx)

    def test_raises_when_no_valid_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            nc = Path(tmp) / "grid.nc"
            x, y = _make_netcdf(nc, nx=5, ny=5)

            ny, nx = 5, 5
            los_e = np.full((ny, nx), np.nan, dtype=np.float32)
            los_n = np.full((ny, nx), np.nan, dtype=np.float32)
            los_u = np.full((ny, nx), np.nan, dtype=np.float32)

            gdf = _make_gnss_gdf(
                x_coords=np.array([x[2]], dtype=np.float32),
                y_coords=np.array([y[2]], dtype=np.float32),
                deast=np.array([0.01]),
                dnorth=np.array([0.005]),
                dup=np.array([0.002]),
            )

            with pytest.raises(ValueError, match="No valid GNSS LOS samples"):
                project_to_los(los_e, los_n, los_u, nc, gdf, method="griddata")

    def test_out_of_swath_fill_is_ignored(self):
        """Regression: `(0, 0, 1)` outside the swath was used as a real look.

        Stations there entered the GNSS field with their vertical rate only,
        biasing it near the swath edge (NYC frame F08622).
        """
        from venti.gnss.los import project_uncertainty_to_los

        with tempfile.TemporaryDirectory() as tmp:
            nc = Path(tmp) / "grid.nc"
            x, y = _make_netcdf(nc, nx=12, ny=8)
            look = (-0.6, -0.1, 0.78)
            los = [np.full((8, 12), c, dtype=np.float32) for c in look]
            los[0][:, 8:], los[1][:, 8:], los[2][:, 8:] = 0.0, 0.0, 1.0  # no swath
            cols, rows = [0, 4, 7, 9, 11], [0, 3, 7]  # cols 9 and 11: no swath
            gdf = _make_gnss_gdf(
                x_coords=np.array([x[c] for c in cols for _ in rows]),
                y_coords=np.array([y[r] for _ in cols for r in rows]),
                deast=np.full(15, -15.0),
                dnorth=np.full(15, 4.0),
                dup=np.full(15, -1.0),
            ).assign(dsigma_e=0.5, dsigma_n=0.5, dsigma_u=2.0)

            inside = gdf[gdf.geometry.x < x[8]]
            for project in (project_to_los, project_uncertainty_to_los):
                np.testing.assert_array_equal(
                    project(*los, nc, gdf, method="rbf"),
                    project(*los, nc, inside, method="rbf"),
                )


@pytest.mark.skipif(not HAS_GEOPANDAS, reason="geopandas not installed")
class TestGNSSReference:
    """Tests for GNSSReference dataclass (no network)."""

    def test_raises_if_utm_epsg_not_set(self):
        from venti.gnss.reference import GNSSReference

        gnss = GNSSReference(
            bounds=(3_800_000, 3_900_000, 400_000, 500_000),
            output_dir=Path("/tmp"),
        )
        with pytest.raises(ValueError, match="utm_epsg"):
            gnss.download_stations()

    def test_default_reference_frame(self):
        from venti.gnss.reference import GNSSReference

        gnss = GNSSReference(
            bounds=(3_800_000, 3_900_000, 400_000, 500_000),
            output_dir=Path("/tmp"),
        )
        assert gnss.reference_frame == "IGS20"

    def test_station_files_empty_before_download(self):
        from venti.gnss.reference import GNSSReference

        gnss = GNSSReference(
            bounds=(3_800_000, 3_900_000, 400_000, 500_000),
            output_dir=Path("/tmp"),
            utm_epsg=32611,
        )
        assert gnss.station_files == []
        assert gnss.station_gdf is None


def _gnss_reference_from_files(tmp_path, grid_type, rates):
    """A GNSSReference over synthetic station files (no download)."""
    from venti.gnss.reference import GNSSReference

    tmp_path.mkdir(parents=True, exist_ok=True)
    years = 2016.0 + np.arange(0, 6 * 365) / 365.25
    files, points = [], []
    for sid, (x, y, ve, vn, vu) in enumerate(rates, start=1):
        path = tmp_path / f"{sid:06d}_IGS20_{grid_type}.tenv8"
        t = years - years[0]
        _write_tenv8(path, years, ve * t, vn * t, vu * t, sigma=0.5)
        files.append(path)
        points.append(Point(x, y))
    gnss = GNSSReference(bounds=(0, 1, 0, 1), output_dir=tmp_path, grid_type=grid_type)
    gnss.station_files = files
    gnss.station_gdf = gpd.GeoDataFrame(index=range(1, len(rates) + 1), geometry=points)
    return gnss


def test_constant_and_variable_grids_give_the_same_los_field(tmp_path):
    """Same linear motion, read as rates or as epoch positions -> same field."""
    from venti.gnss import compute_gnss_los, compute_gnss_los_std

    rates = [  # x, y (m), then ve, vn, vu (mm/yr)
        (100.0, 500.0, -14.0, -3.0, -1.0),
        (700.0, 450.0, -13.0, -3.5, 0.5),
        (400.0, 100.0, -14.5, -2.5, -2.0),
        (650.0, 150.0, -13.5, -2.8, 0.0),
    ]
    x, y = np.arange(0.0, 800, 20), np.arange(600.0, 0, -20)
    los = [np.full((y.size, x.size), c) for c in (-0.6, -0.1, 0.78)]
    ref_date, sec_date = 2017.0, 2019.5

    fields = {}
    for grid_type in ("constant", "variable"):
        gnss = _gnss_reference_from_files(tmp_path / grid_type, grid_type, rates)
        fields[grid_type] = compute_gnss_los(
            gnss, *los, (x, y), ref_date=ref_date, sec_date=sec_date
        )
        std = compute_gnss_los_std(
            gnss, *los, (x, y), ref_date=ref_date, sec_date=sec_date
        )
        assert np.isfinite(std).all()
        assert (std > 0).all()

    # Positions are daily, so dates round to the nearest day (< 0.1 mm here).
    np.testing.assert_allclose(fields["variable"], fields["constant"], atol=0.1)
    # sec - ref of a westward-moving network seen with los_e < 0: positive.
    assert np.nanmedian(fields["constant"]) > 0
