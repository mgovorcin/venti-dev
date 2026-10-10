# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Unit tests for venti.workflow.stage_frame_data.

Covers pure logic and orchestration — no network or file I/O is exercised.
All external calls (downloads, DEM/LOS generation, tropo processing, GNSS
downloads) are patched with mocks.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

try:
    import geopandas as gpd
    from shapely.geometry import Point

    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False

from venti.workflow.stage_frame_data import _find_disp_file, download_gnss_data

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_DISP_NAME = "OPERA_L3_DISP-S1_IW_F08887_VV_20160705T002809Z_20160616T002812Z_v1.0.nc"


def _make_disp_file(directory: Path, name: str = _DISP_NAME) -> Path:
    """Write an empty file with a DISP-S1-style name."""
    path = directory / name
    path.touch()
    return path


def _make_station_gdf(station_ids: list[int], epsg: int = 32611) -> gpd.GeoDataFrame:
    """Build a minimal station GeoDataFrame for a given list of IDs."""
    return gpd.GeoDataFrame(
        {"lon": [-118.0] * len(station_ids), "lat": [35.0] * len(station_ids)},
        index=station_ids,
        geometry=[Point(400_000 + i * 5_000, 3_850_000) for i in station_ids],
        crs=f"EPSG:{epsg}",
    )


def _write_tenv8(path: Path, n: int = 30) -> None:
    """Write a minimal .tenv8 file with a linear trend."""
    years = np.linspace(2014.0, 2024.0, n)
    rows = np.column_stack(
        [
            years,
            years * 0.010,
            years * -0.005,
            years * 0.002,
            np.full(n, 1e-3),
            np.full(n, 1e-3),
            np.full(n, 1e-3),
        ]
    )
    np.savetxt(path, rows, fmt="%.6f")


# ---------------------------------------------------------------------------
# _find_disp_file
# ---------------------------------------------------------------------------


class TestFindDispFile:
    """Tests for _find_disp_file."""

    def test_returns_single_matching_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            disp_dir = Path(tmp)
            expected = _make_disp_file(disp_dir)
            result = _find_disp_file(disp_dir, frame_id=8887)
            assert result == expected

    def test_raises_when_no_files(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            pytest.raises(FileNotFoundError, match="No DISP-S1 file found"),
        ):
            _find_disp_file(Path(tmp), frame_id=8887)

    def test_raises_when_multiple_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            disp_dir = Path(tmp)
            _make_disp_file(
                disp_dir,
                "OPERA_L3_DISP-S1_IW_F08887_VV_20160101T000000Z_20160601T000000Z_v1.0.nc",
            )
            _make_disp_file(
                disp_dir,
                "OPERA_L3_DISP-S1_IW_F08887_VV_20160101T000000Z_20161201T000000Z_v1.0.nc",
            )
            with pytest.raises(ValueError, match="Multiple DISP-S1 files found"):
                _find_disp_file(disp_dir, frame_id=8887)

    def test_ignores_different_frame_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            disp_dir = Path(tmp)
            _make_disp_file(disp_dir)  # frame 8887
            with pytest.raises(FileNotFoundError):
                _find_disp_file(disp_dir, frame_id=9999)


# ---------------------------------------------------------------------------
# download_gnss_data
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_GEOPANDAS, reason="geopandas not installed")
class TestDownloadGnssData:
    """Tests for download_gnss_data — all network calls are mocked."""

    _FRAME_ID = 8887
    _EPSG = 32611

    def _bbox(self):
        """Return a minimal bbox-like object."""
        bbox = MagicMock()
        bbox.bottom = 3_800_000.0
        bbox.top = 3_900_000.0
        bbox.left = 400_000.0
        bbox.right = 500_000.0
        return bbox

    def test_velocities_parquet_written(self):
        station_gdf = _make_station_gdf([1, 2])

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "gnss"
            stations_dir = output_dir / "stations"

            # Pre-write station files so calculate_station_velocity has real input
            output_dir.mkdir(parents=True)
            stations_dir.mkdir()
            for sid in [1, 2]:
                _write_tenv8(stations_dir / f"{sid:06d}_IGS20.tenv8")

            with (
                patch(
                    "venti.workflow.stage_frame_data.get_frame_bbox",
                    return_value=(self._EPSG, self._bbox()),
                ),
                patch(
                    "venti.workflow.stage_frame_data.download_grid_lookup",
                    return_value=output_dir / "grid_latlon_lookup.txt",
                ),
                patch(
                    "venti.workflow.stage_frame_data.find_stations_in_bounds",
                    return_value=station_gdf,
                ),
                patch(
                    "venti.workflow.stage_frame_data.download_station",
                    side_effect=lambda sid, d, **_: d / f"{sid:06d}_IGS20.tenv8",
                ),
            ):
                vel_path = download_gnss_data(self._FRAME_ID, output_dir, num_workers=1)

            assert vel_path.exists()
            df = gpd.read_parquet(vel_path)
            assert len(df) == 2
            assert {"ve_mmyr", "vn_mmyr", "vu_mmyr"}.issubset(df.columns)

    def test_raises_when_no_stations_in_bounds(self):
        empty_gdf = _make_station_gdf([])

        with (
            tempfile.TemporaryDirectory() as tmp,
            patch(
                "venti.workflow.stage_frame_data.get_frame_bbox",
                return_value=(self._EPSG, self._bbox()),
            ),
            patch(
                "venti.workflow.stage_frame_data.download_grid_lookup",
                return_value=Path(tmp) / "grid_latlon_lookup.txt",
            ),
            patch(
                "venti.workflow.stage_frame_data.find_stations_in_bounds",
                return_value=empty_gdf,
            ),
            pytest.raises(ValueError, match="No UNR stations found"),
        ):
            download_gnss_data(self._FRAME_ID, Path(tmp) / "gnss")

    def test_failed_station_download_is_skipped(self):
        """A RuntimeError on one station should not abort the whole run."""
        station_gdf = _make_station_gdf([1, 2])

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "gnss"
            stations_dir = output_dir / "stations"
            output_dir.mkdir(parents=True)
            stations_dir.mkdir()
            _write_tenv8(stations_dir / "000002_IGS20.tenv8")

            def _download_side_effect(sid, d, **_):
                if sid == 1:
                    msg = "HTTP 404"
                    raise RuntimeError(msg)
                return d / f"{sid:06d}_IGS20.tenv8"

            with (
                patch(
                    "venti.workflow.stage_frame_data.get_frame_bbox",
                    return_value=(self._EPSG, self._bbox()),
                ),
                patch(
                    "venti.workflow.stage_frame_data.download_grid_lookup",
                    return_value=output_dir / "grid_latlon_lookup.txt",
                ),
                patch(
                    "venti.workflow.stage_frame_data.find_stations_in_bounds",
                    return_value=station_gdf,
                ),
                patch(
                    "venti.workflow.stage_frame_data.download_station",
                    side_effect=_download_side_effect,
                ),
            ):
                vel_path = download_gnss_data(self._FRAME_ID, output_dir, num_workers=1)

            df = gpd.read_parquet(vel_path)
            assert len(df) == 1

    def test_raises_when_all_velocity_estimations_fail(self):
        station_gdf = _make_station_gdf([1])

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "gnss"
            output_dir.mkdir(parents=True)
            (output_dir / "stations").mkdir()

            with (
                patch(
                    "venti.workflow.stage_frame_data.get_frame_bbox",
                    return_value=(self._EPSG, self._bbox()),
                ),
                patch(
                    "venti.workflow.stage_frame_data.download_grid_lookup",
                    return_value=output_dir / "grid_latlon_lookup.txt",
                ),
                patch(
                    "venti.workflow.stage_frame_data.find_stations_in_bounds",
                    return_value=station_gdf,
                ),
                patch(
                    "venti.workflow.stage_frame_data.download_station",
                    side_effect=lambda sid, d, **_: d / f"{sid:06d}_IGS20.tenv8",
                ),
                patch(
                    "venti.workflow.stage_frame_data.calculate_station_velocity",
                    side_effect=ValueError("bad file"),
                ),
                pytest.raises(RuntimeError, match="Velocity estimation failed"),
            ):
                download_gnss_data(self._FRAME_ID, output_dir, num_workers=1)


# ---------------------------------------------------------------------------
# stage_frame orchestration
# ---------------------------------------------------------------------------


def _make_staging_mocks(tmp: Path) -> dict:
    """Build a dict of mock staging functions and their expected return paths."""
    disp_dir = tmp / "disp_s1"
    disp_dir.mkdir(parents=True)
    disp_file = _make_disp_file(disp_dir)

    dem_path = tmp / "dem" / "dem_frame_8887.tif"
    los_path = tmp / "los" / "los_enu_frame_8887.tif"
    inc_path = tmp / "los" / "incidence_angle_frame_8887.tif"

    return {
        "disp_file": disp_file,
        "dem_path": dem_path,
        "los_path": los_path,
        "inc_path": inc_path,
        "download_frame_products": MagicMock(),
        "generate_frame_dem": MagicMock(return_value=dem_path),
        "generate_los_enu_raster": MagicMock(return_value=los_path),
        "generate_incidence_angle_raster": MagicMock(return_value=inc_path),
        "process_tropo_from_file": MagicMock(),
        "download_gnss_data_fn": MagicMock(
            return_value=tmp / "gnss" / "velocities.parquet"
        ),
    }


def _inject_staging_modules(mocks: dict) -> dict[str, ModuleType]:
    """Return fake modules to inject into sys.modules for staging CLI imports."""
    dem_mod = ModuleType("dem_cli")
    dem_mod.generate_frame_dem = mocks["generate_frame_dem"]  # type: ignore[attr-defined]

    disp_mod = ModuleType("disp_cli")
    disp_mod.download_frame_products = mocks["download_frame_products"]  # type: ignore[attr-defined]

    los_mod = ModuleType("los_cli")
    los_mod.generate_los_enu_raster = mocks["generate_los_enu_raster"]  # type: ignore[attr-defined]
    los_mod.generate_incidence_angle_raster = mocks["generate_incidence_angle_raster"]  # type: ignore[attr-defined]

    tropo_mod = ModuleType("tropo_cli")
    tropo_mod.process_tropo_from_file = mocks["process_tropo_from_file"]  # type: ignore[attr-defined]

    utils_mod = ModuleType("utils")
    utils_mod.parse_date = MagicMock(return_value=MagicMock())  # type: ignore[attr-defined]

    return {
        "dem_cli": dem_mod,
        "disp_cli": disp_mod,
        "los_cli": los_mod,
        "tropo_cli": tropo_mod,
        "utils": utils_mod,
    }


class TestStageFrame:
    """Tests for stage_frame orchestration — all I/O is mocked."""

    def _run(
        self,
        tmp: Path,
        skip_tropo: bool = False,
        skip_gnss: bool = False,
    ) -> dict:
        mocks = _make_staging_mocks(tmp)
        fake_modules = _inject_staging_modules(mocks)

        with (
            patch.dict("sys.modules", fake_modules),
            patch(
                "venti.workflow.stage_frame_data.download_gnss_data",
                mocks["download_gnss_data_fn"],
            ),
        ):
            from venti.workflow.stage_frame_data import stage_frame

            stage_frame(
                frame_id=8887,
                date="2016-06-15",
                output_dir=tmp,
                skip_tropo=skip_tropo,
                skip_gnss=skip_gnss,
            )

        return mocks

    def test_all_steps_called_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            mocks = self._run(Path(tmp))

        mocks["download_frame_products"].assert_called_once()
        mocks["generate_frame_dem"].assert_called_once()
        mocks["generate_los_enu_raster"].assert_called_once()
        mocks["generate_incidence_angle_raster"].assert_called_once()
        mocks["process_tropo_from_file"].assert_called_once()
        mocks["download_gnss_data_fn"].assert_called_once()

    def test_skip_tropo_omits_tropo_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            mocks = self._run(Path(tmp), skip_tropo=True)

        mocks["process_tropo_from_file"].assert_not_called()
        mocks["download_gnss_data_fn"].assert_called_once()

    def test_skip_gnss_omits_gnss_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            mocks = self._run(Path(tmp), skip_gnss=True)

        mocks["process_tropo_from_file"].assert_called_once()
        mocks["download_gnss_data_fn"].assert_not_called()

    def test_skip_both_runs_core_steps_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            mocks = self._run(Path(tmp), skip_tropo=True, skip_gnss=True)

        mocks["download_frame_products"].assert_called_once()
        mocks["generate_frame_dem"].assert_called_once()
        mocks["generate_los_enu_raster"].assert_called_once()
        mocks["generate_incidence_angle_raster"].assert_called_once()
        mocks["process_tropo_from_file"].assert_not_called()
        mocks["download_gnss_data_fn"].assert_not_called()

    def test_dem_uses_correct_output_subdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            mocks = self._run(tmp_path, skip_tropo=True, skip_gnss=True)

        _, kwargs = mocks["generate_frame_dem"].call_args
        assert kwargs["output_dir"] == tmp_path / "dem"

    def test_tropo_receives_dem_and_incidence_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            mocks = self._run(tmp_path, skip_gnss=True)

        _, kwargs = mocks["process_tropo_from_file"].call_args
        assert kwargs["dem_path"] == mocks["dem_path"]
        assert kwargs["incidence_angle_path"] == mocks["inc_path"]
