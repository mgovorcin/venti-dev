# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tests for venti.models.load_itrf.

`convert_to_euler_poles` imported from a module path that does not exist
(`plate_motion.euler_pole`), so it raised ImportError on first use. These
tests pin the import and check the conversion against the function it wraps.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from venti.models import load_itrf
from venti.models.plate_motion import rotation_rate_to_euler_pole

# ITRF2020-PMM North America rotation rate, deg/Myr (Altamimi et al. 2023,
# rounded): a pole near 88 W, 5 S, rotating at ~0.70 deg/Myr.
NOAM = {"omega_x": 0.024, "omega_y": -0.694, "omega_z": -0.063}


def test_convert_to_euler_poles_matches_plate_motion():
    df = pd.DataFrame({"plate": ["NOAM"], "name": ["North America"], **{k: [v] for k, v in NOAM.items()}})
    out = load_itrf.convert_to_euler_poles(df)

    assert list(out.columns) == ["plate", "euler_longitude", "euler_latitude", "angular_velocity", "name"]
    assert out.loc[0, "plate"] == "NOAM"
    assert out.loc[0, "name"] == "North America"

    to_rad_per_yr = np.deg2rad(1.0) / 1e6
    lon, lat, omega = rotation_rate_to_euler_pole(
        NOAM["omega_x"] * to_rad_per_yr,
        NOAM["omega_y"] * to_rad_per_yr,
        NOAM["omega_z"] * to_rad_per_yr,
    )
    assert out.loc[0, "euler_longitude"] == pytest.approx(lon)
    assert out.loc[0, "euler_latitude"] == pytest.approx(lat)
    assert out.loc[0, "angular_velocity"] == pytest.approx(omega)

    # Physical sanity, independent of the wrapped function.
    assert out.loc[0, "euler_longitude"] == pytest.approx(-88.0, abs=1.0)
    assert out.loc[0, "euler_latitude"] == pytest.approx(-5.2, abs=0.5)
    assert out.loc[0, "angular_velocity"] == pytest.approx(0.697, abs=0.005)


def test_convert_to_euler_poles_without_name_column():
    df = pd.DataFrame({"plate": ["NOAM"], **{k: [v] for k, v in NOAM.items()}})
    out = load_itrf.convert_to_euler_poles(df)
    assert "name" not in out.columns
    assert len(out) == 1


def test_convert_to_euler_poles_requires_rotation_columns():
    with pytest.raises(KeyError):
        load_itrf.convert_to_euler_poles(pd.DataFrame({"plate": ["NOAM"]}))
