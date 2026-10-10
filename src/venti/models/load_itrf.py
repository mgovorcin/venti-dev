# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# Get the directory where this module is located
MODULE_DIR = Path(__file__).parent
ITRF14_DATA = MODULE_DIR / "data/itrf2014_sites.csv"
ITRF20_DATA = MODULE_DIR / "data/itrf2020_sites.csv"


def json_to_dataframe(json_data: dict) -> pd.DataFrame:
    """Convert JSON-format ITRF Plate Motion Model data to a pandas DataFrame.

    Parameters
    ----------
    json_data : dict
        Dictionary containing ITRF PMM data from `load_itrf_pmm()`.

    Returns
    -------
    pandas.DataFrame
            DataFrame with columns: 'plate', 'name', 'omega_x', 'omega_y', 'omega_z'.

    """
    rows = []
    for plate_code, plate_data in json_data["plates"].items():
        rows.append(
            {
                "plate": plate_code,
                "name": plate_data["name"],
                "omega_x": plate_data["omega_x"],
                "omega_y": plate_data["omega_y"],
                "omega_z": plate_data["omega_z"],
            }
        )

    return pd.DataFrame(rows)


def load_itrf_pmm(filepath: str | Path | None = None, date: int = 2020) -> dict:
    """Load ITRF Plate Motion Model (PMM) data from JSON format.

    Parameters
    ----------
    filepath : str or Path, optional
        Path to the JSON file containing ITRF PMM data. If None, the default
        file `data/itrf{date}_pmm.json` in the module directory is used.
    date : int, optional
        ITRF reference frame year (e.g., 2014, 2020). Default is 2020.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing plate rotation parameters.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If the JSON file is invalid or malformed.

    """
    if filepath is None:
        filepath = MODULE_DIR / f"data/itrf{date}_pmm.json"

    try:
        with open(filepath) as f:
            data = json.load(f)

        # Basic validation
        if "metadata" not in data or "plates" not in data:
            msg = f"Invalid ITRF JSON format in {filepath}"
            raise ValueError(msg)

        return json_to_dataframe(data)

    except FileNotFoundError as e:
        msg = f"ITRF data file not found: {filepath}"
        raise FileNotFoundError(msg) from e
    except json.JSONDecodeError as e:
        msg = f"Invalid JSON in {filepath}: {e}"
        raise ValueError(msg) from e


def convert_to_euler_poles(df: pd.DataFrame) -> pd.DataFrame:
    """Convert ITRF rotation rates to Euler pole parameters.

    Converts plate rotation rates (`omega_x`, `omega_y`, `omega_z`) from
    degrees per million years (deg/Myr) to Euler pole representation
    (longitude, latitude, angular velocity).

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns 'omega_x', 'omega_y', 'omega_z' (deg/Myr),
        and optionally 'plate' and 'name'.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns 'plate', 'euler_longitude', 'euler_latitude',
        'angular_velocity', and optionally 'name'.

    """
    from .plate_motion import rotation_rate_to_euler_pole

    results = []

    for _, row in df.iterrows():
        # Convert from deg/Myr to rad/year
        wx = np.radians(row["omega_x"]) / 1e6  # deg/Myr → rad/year
        wy = np.radians(row["omega_y"]) / 1e6
        wz = np.radians(row["omega_z"]) / 1e6

        # Convert to Euler pole
        euler_lon, euler_lat, omega = rotation_rate_to_euler_pole(wx, wy, wz)

        result = {
            "plate": row["plate"],
            "euler_longitude": euler_lon,
            "euler_latitude": euler_lat,
            "angular_velocity": omega,  # deg/Myr
        }

        # Include plate name if available
        if "name" in row:
            result["name"] = row["name"]

        results.append(result)

    return pd.DataFrame(results)


def get_plate_data(plate_code: str, date: int = 2020) -> dict:
    """Retrieve rotation data for a specific tectonic plate.

    Parameters
    ----------
    plate_code : str
        Plate abbreviation (e.g., 'PCFC', 'NOAM').
    date : int, optional
        ITRF reference frame year. Default is 2020.

    Returns
    -------
    dict
        Dictionary containing plate rotation parameters including
        'plate', 'euler_longitude', 'euler_latitude', 'angular_velocity',
        and optionally 'name'.

    Raises
    ------
    ValueError
        If the specified plate is not found in the dataset.

    """
    df = load_itrf_pmm(date=date)

    plate_row = df[df["plate"] == plate_code.upper()]
    if len(plate_row) == 0:
        available_plates = df["plate"].tolist()
        msg = f"Plate '{plate_code}' not found. Available: {available_plates}"
        raise ValueError(msg)

    return plate_row.iloc[0].to_dict()
