# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Plate motion is geepers' (ADR-0010/0020); venti.models only re-exports it."""

from __future__ import annotations

import numpy as np
import pytest

import venti.models as models


def test_models_reexports_geepers_plate_motion():
    from geepers import euler

    assert models.plate_velocity_enu is euler.plate_velocity_enu
    assert models.plate_pole is euler.plate_pole
    assert models.PLATE_CODES["NA"] == "NOAM"


def test_removed_modules_are_gone():
    with pytest.raises(ImportError):
        import venti.models.plate_motion
    with pytest.raises(ImportError):
        import venti.models.load_itrf  # noqa: F401


def test_get_model_rates_uses_geepers(monkeypatch):
    """The PMM raster helper produces mm/yr grids from the geepers pole."""
    gpd = pytest.importorskip("geopandas")
    import importlib.util
    from pathlib import Path

    from shapely.geometry import box

    spec = importlib.util.spec_from_file_location(
        "get_model_rates",
        Path(__file__).resolve().parents[1] / "scripts/get_model_rates.py",
    )
    gmr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gmr)
    frame = gpd.GeoDataFrame(geometry=[box(-96.0, 29.0, -95.0, 30.0)], crs="EPSG:4326")
    ve, vn, attrs = gmr.get_frame_pmm(
        frame, plate="NA", date=2020, grid_posting=10000  # frame in EPSG:4326, posting in metres
    )
    assert ve.shape == vn.shape == (attrs["height"], attrs["width"])
    assert attrs["units"] == "mm/year"
    assert attrs["model"] == "ITRF2020"
    # Houston area, NA in ITRF2020: ~12-13 mm/yr west, a few mm/yr south
    assert np.all((ve > -16) & (ve < -9))
    assert np.all(vn < 0)
