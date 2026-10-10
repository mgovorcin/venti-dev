# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Euler Pole Calculation for Tectonic Plate Motion.

This module provides functions to calculate Euler poles from GPS velocity data,
handling coordinate transformations and uncertainty estimation.
"""

from __future__ import annotations

import numpy as np
from pyproj import Transformer

# Constants
EARTH_RADIUS_KM = 6371.0
DEG_TO_RAD = np.pi / 180.0
RAD_TO_DEG = 180.0 / np.pi
MYR_TO_YEAR = 1e6
MAS_TO_RAD = np.pi / 648000000


# Initialize coordinate transformer (WGS84 to geocentric cartesian)
LONLAT_TO_XYZ_TRANSFORMER = Transformer.from_crs(
    "EPSG:4326",
    {"proj": "geocent", "ellps": "GRS80", "datum": "WGS84"},
    always_xy=True,
)


def get_conversion_matrix(longitude: float, latitude: float) -> np.ndarray:
    """Compute the conversion matrix from geocentric XYZ to local ENU.

    Parameters
    ----------
    longitude : float
        Longitude of the site in degrees.
    latitude : float
        Latitude of the site in degrees.

    Returns
    -------
    np.ndarray
        3x3 rotation matrix converting XYZ → ENU coordinates.

    """
    lon_rad = np.radians(longitude)
    lat_rad = np.radians(latitude)

    cos_lon, sin_lon = np.cos(lon_rad), np.sin(lon_rad)
    cos_lat, sin_lat = np.cos(lat_rad), np.sin(lat_rad)

    R = np.array(
        [
            [-sin_lon, cos_lon, 0],
            [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
            [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
        ],
        dtype=np.float64,
    )

    return R


def get_local_frame(
    longitude: float, latitude: float, height: float = 0.0
) -> np.ndarray:
    """Compute the local observation frame (RAi) for a site.

    Parameters
    ----------
    longitude : float
        Longitude in degrees.
    latitude : float
        Latitude in degrees.
    height : float, optional
        Height above ellipsoid in meters, by default 0.0

    Returns
    -------
    np.ndarray
        3x3 local observation matrix relating rotation rate to site velocity.

    """
    # Convert to geocentric cartesian coordinates
    x, y, z = LONLAT_TO_XYZ_TRANSFORMER.transform(
        longitude, latitude, height, radians=False
    )

    # Get conversion matrix XYZ to ENU
    R = get_conversion_matrix(longitude, latitude)

    # Observation equation matrix (skew-symmetric for cross product)
    # ω × r = [ωy*z - ωz*y, ωz*x - ωx*z, ωx*y - ωy*x]  # noqa: RUF003
    Ai = np.array([[0, z, -y], [-z, 0, x], [y, -x, 0]], dtype=np.float64)

    # Apply rotation matrix
    RAi = np.dot(R, Ai)

    return RAi


def rotation_rate_to_euler_pole(
    wx: float, wy: float, wz: float
) -> tuple[float, float, float]:
    """Convert rotation rate vector to Euler pole parameters.

    Parameters
    ----------
    wx : float
        Rotation rate component along x-axis in rad/year.
    wy : float
        Rotation rate component along y-axis in rad/year.
    wz : float
        Rotation rate component along z-axis in rad/year.

    Returns
    -------
    tuple of float
        (euler_longitude, euler_latitude, omega)
        - euler_longitude : degrees
        - euler_latitude : degrees
        - omega : angular velocity in deg/Myr

    """
    W_magnitude = np.sqrt(wx**2 + wy**2 + wz**2)

    # Angular velocity in deg/Myr from rad/year
    # rad/year × (180°/π rad) × (1e6 year/Myr) = deg/Myr  # noqa: RUF003
    omega = np.degrees(W_magnitude) * MYR_TO_YEAR

    # Euler pole latitude (colatitude from z-axis)
    euler_latitude = 90.0 - np.degrees(np.arccos(wz / W_magnitude))

    # Euler pole longitude
    euler_longitude = np.degrees(np.arctan2(wy, wx))

    return euler_longitude, euler_latitude, omega


def euler_pole_to_rotation_rate(
    euler_longitude: float, euler_latitude: float, omega: float
) -> tuple[float, float, float]:
    """Convert Euler pole parameters to rotation rate vector.

    Parameters
    ----------
    euler_longitude : float
        Euler pole longitude in degrees.
    euler_latitude : float
        Euler pole latitude in degrees.
    omega : float
        Angular velocity in deg/Myr.

    Returns
    -------
    tuple of float
        Rotation rate components (wx, wy, wz) in rad/year.

    """
    # Convert angular velocity from deg/Myr to rad/year
    # deg/Myr × (π rad/180°) × (1 Myr/1e6 year) = rad/year  # noqa: RUF003
    omega_rad_per_year = np.radians(omega) / MYR_TO_YEAR

    cos_lon = np.cos(np.radians(euler_longitude))
    sin_lon = np.sin(np.radians(euler_longitude))
    cos_lat = np.cos(np.radians(euler_latitude))
    sin_lat = np.sin(np.radians(euler_latitude))

    wx = cos_lat * cos_lon * omega_rad_per_year
    wy = cos_lat * sin_lon * omega_rad_per_year
    wz = sin_lat * omega_rad_per_year

    return wx, wy, wz


def calculate_euler_pole(
    longitude: np.ndarray,
    latitude: np.ndarray,
    velocity_east: np.ndarray,
    velocity_north: np.ndarray,
    sigma_east: np.ndarray,
    sigma_north: np.ndarray,
    heights: np.ndarray | None = None,
    correlation_coefficient: float | None = 0,
) -> tuple[float, float, float, dict]:
    """Calculate the Euler pole parameters from GPS velocity observations.

    Parameters
    ----------
    longitude : ndarray
        Site longitudes in degrees. Shape: (n_sites,)
    latitude : ndarray
        Site latitudes in degrees. Shape: (n_sites,)
    velocity_east : ndarray
        East velocity components in mm/year. Shape: (n_sites,)
    velocity_north : ndarray
        North velocity components in mm/year. Shape: (n_sites,)
    sigma_east : ndarray
        Uncertainties of east velocity components in mm/year. Shape: (n_sites,)
    sigma_north : ndarray
        Uncertainties of north velocity components in mm/year. Shape: (n_sites,)
    heights : ndarray, optional
        Heights above the ellipsoid in meters. If None, all heights are assumed 0.
        Shape: (n_sites,)
    correlation_coefficient : float, optional
        Correlation coefficient between east and north velocity components. Default
        is 0.

    Returns
    -------
    euler_longitude : float
        Longitude of the Euler pole in degrees.
    euler_latitude : float
        Latitude of the Euler pole in degrees.
    omega : float
        Angular velocity of the Euler pole in deg/Myr.
    statistics : dict
        Dictionary containing fit statistics and quality metrics:
        - 'rms' : float
            Root mean square residual in mm/year
        - 'wrms' : float
            Weighted RMS residual in mm/year
        - 'chi_squared' : float
            Chi-squared value
        - 'reduced_chi_squared' : float
            Reduced chi-squared value
        - 'degrees_of_freedom' : int
            Number of degrees of freedom
        - 'parameter_covariance' : ndarray
            3x3 covariance matrix of rotation rates in (rad/year)^2
        - 'rotation_rates_rad_per_year' : list of float
            [wx, wy, wz] rotation rates in rad/year

    Raises
    ------
    ValueError
        If the covariance matrix or normal matrix is singular, degrees of freedom are
        insufficient, or covariance matrix is not positive definite.

    Notes
    -----
    - The method applies weighted least squares to estimate rotation parameters from
    observed horizontal velocities at multiple GPS sites.
    - Rotation rates (wx, wy, wz) are internally converted from deg/Myr to rad/year.
    - The returned Euler pole parameters are in geodetic coordinates (lon/lat) and
    deg/Myr.
    - The function propagates observational uncertainties and optionally accounts
    for correlation between east and north components.

    """
    n_sites = len(longitude)

    if heights is None:
        heights = np.zeros(n_sites)

    # Initialize design matrix and observation vector
    A = np.zeros((2 * n_sites, 3), dtype=np.float64)
    b = np.zeros((2 * n_sites, 1), dtype=np.float64)
    cov = np.zeros((2 * n_sites, 2 * n_sites), dtype=np.float64)

    # Build system matrices
    for i, (lon, lat, hgt, ve, vn, se, sn) in enumerate(
        zip(
            longitude,
            latitude,
            heights,
            velocity_east,
            velocity_north,
            sigma_east,
            sigma_north,
            strict=False,
        )
    ):
        # Get local frame matrix
        local_frame = get_local_frame(lon, lat, hgt)

        # Correlation coefficient should be dimensionless between -1 and 1
        # Covariance = σe × σn × correlation_coefficient  # noqa: RUF003
        cross_correlation = se * sn * correlation_coefficient

        # Fill design matrix
        A[2 * i, :] = local_frame[0, :]  # East component
        A[2 * i + 1, :] = local_frame[1, :]  # North component

        # Fill observation vector
        b[2 * i, 0] = ve
        b[2 * i + 1, 0] = vn

        # Fill covariance matrix
        cov[2 * i, 2 * i] = se**2
        cov[2 * i + 1, 2 * i + 1] = sn**2
        cov[2 * i + 1, 2 * i] = cross_correlation
        cov[2 * i, 2 * i + 1] = cross_correlation

    # Weighted least squares solution
    eigenvals = np.linalg.eigvals(cov)
    if np.any(eigenvals <= 0):
        msg = "Covariance matrix is not positive definite"
        raise ValueError(msg)

    # Solve linear system

    try:
        P = np.linalg.inv(cov)
    except np.linalg.LinAlgError as err:
        msg = "Covariance matrix is singular."
        raise ValueError(msg) from err

    ATP = A.T @ P

    N = ATP @ A
    M = ATP @ b

    # Solve for rotation parameters
    try:
        Q = np.linalg.inv(N)
    except np.linalg.LinAlgError as err:
        msg = "Normal matrix is singular."
        raise ValueError(msg) from err

    X = np.dot(Q, M)

    # Calculate model predictions and residuals
    MP = np.dot(A, X)
    residuals = b - MP

    # Statistical analysis
    chi2 = float(residuals.T @ P @ residuals)
    # Degrees of freedom
    dof = 2 * n_sites - 3
    if dof <= 0:
        msg = f"Insufficient degrees of freedom: {dof}. Need at least 3 sites."
        raise ValueError(msg)

    reduced_chi2 = np.sqrt(chi2 / dof)

    # Extract rotation components
    wx = X[0, 0] / 1000  # rad/year
    wy = X[1, 0] / 1000  # rad/year
    wz = X[2, 0] / 1000  # rad/year

    # Convert to Euler pole
    euler_longitude, euler_latitude, omega = rotation_rate_to_euler_pole(wx, wy, wz)

    # Calculate RMS statistics
    re = residuals[::2].flatten()  # residuals east
    rn = residuals[1::2].flatten()  # residuals north

    rms = np.sqrt(np.mean(re**2 + rn**2))

    # Weighted RMS
    try:
        wrms = np.sum((re / sigma_east) ** 2 + (rn / sigma_north) ** 2)
        wrms /= np.sum(1 / sigma_east**2 + 1 / sigma_north**2)
        wrms = np.sqrt(wrms)
    except ZeroDivisionError:
        wrms = np.nan

    # Compile statistics
    statistics = {
        "rms": rms,
        "wrms": wrms,
        "chi_squared": chi2,
        "reduced_chi_squared": reduced_chi2,
        "degrees_of_freedom": dof,
        "parameter_covariance": Q / (1000.0**2),  # (mm/year per km)² to (rad/year)²
        "rotation_rates_rad_per_year": [wx, wy, wz],
    }

    return euler_longitude, euler_latitude, omega, statistics


def get_euler_pole_uncertainty(
    euler_longitude: float,
    euler_latitude: float,
    omega: float,
    parameter_covariance: np.ndarray,
) -> tuple[float, float, float, float]:
    """Calculate the uncertainty of an Euler pole and its angular velocity.

    Parameters
    ----------
    euler_longitude : float
        Euler pole longitude in degrees.
    euler_latitude : float
        Euler pole latitude in degrees.
    omega : float
        Angular velocity in deg/Myr.
    parameter_covariance : ndarray
        3x3 covariance matrix of rotation rate parameters in (rad/year)^2.

    Returns
    -------
    max_sigma : float
        Maximum angular uncertainty of the Euler pole in degrees.
    min_sigma : float
        Minimum angular uncertainty of the Euler pole in degrees.
    azimuth : float
        Orientation of the uncertainty ellipse in degrees.
    sigma_omega : float
        Angular velocity uncertainty in deg/Myr.

    """
    # Get rotation matrix for Euler pole location
    rotation_matrix = get_conversion_matrix(euler_longitude, euler_latitude)

    # Convert Euler pole to rotation rates
    wx, wy, wz = euler_pole_to_rotation_rate(euler_longitude, euler_latitude, omega)
    rotation_magnitude = np.sqrt(wx**2 + wy**2 + wz**2)

    # Transform covariance to local ENU frame
    local_covariance = rotation_matrix @ parameter_covariance @ rotation_matrix.T

    # Extract horizontal components
    cov_11 = local_covariance[0, 0]  # East-East
    cov_12 = local_covariance[0, 1]  # East-North
    cov_22 = local_covariance[1, 1]  # North-North

    # Calculate uncertainty ellipse parameters
    discriminant = np.sqrt((cov_11 - cov_22) ** 2 + 4 * cov_12**2)

    max_eigenvalue = 0.5 * (cov_11 + cov_22 + discriminant)
    min_eigenvalue = 0.5 * (cov_11 + cov_22 - discriminant)

    # Convert to angular uncertainties
    max_sigma = np.degrees(np.arctan(np.sqrt(max_eigenvalue) / rotation_magnitude))
    min_sigma = np.degrees(np.arctan(np.sqrt(min_eigenvalue) / rotation_magnitude))

    # Ellipse orientation
    if cov_11 != cov_22:
        azimuth = np.degrees(0.5 * np.arctan(2 * cov_12 / (cov_11 - cov_22)))
    else:
        azimuth = 0.0

    # Angular velocity uncertainty
    sigma_omega = np.degrees(np.sqrt(local_covariance[2, 2])) * MYR_TO_YEAR

    return max_sigma, min_sigma, azimuth, sigma_omega


def model_plate_velocities(
    longitude: float | np.ndarray,
    latitude: float | np.ndarray,
    wx: float,
    wy: float,
    wz: float,
    heights: float | np.ndarray | None = None,
    rotation_covariance: np.ndarray | None = None,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Compute modeled plate velocities at specified locations from rotation rates.

    Parameters
    ----------
    longitude : float or ndarray
        Site longitude(s) in degrees.
    latitude : float or ndarray
        Site latitude(s) in degrees.
    wx, wy, wz : float
        Rotation rate components in deg/Myr.
    heights : float or ndarray, optional
        Site heights above ellipsoid in meters. If a single value is provided,
        it is applied to all sites. Default is None (0 m).
    rotation_covariance : ndarray, optional
        3x3 covariance matrix of rotation rates in (rad/year)^2. If provided,
        velocity uncertainties are calculated.

    Returns
    -------
    ndarray or tuple of ndarray
        If `rotation_covariance` is None:
            Array of modeled velocities [ve, vn] in mm/year (shape: (n_sites, 2)
            or (2,) for single site).
        If `rotation_covariance` is provided:
            Tuple of (modeled_velocities, velocity_uncertainties)
            - modeled_velocities : ndarray of shape (n_sites, 2) in mm/year
            - velocity_uncertainties : ndarray of shape (n_sites, 2) in mm/year

    """
    # Convert inputs to arrays for consistent handling
    longitude = np.atleast_1d(longitude)
    latitude = np.atleast_1d(latitude)

    if len(longitude) != len(latitude):
        msg = "Longitude and latitude arrays must have the same length"
        raise ValueError(msg)

    n_sites = len(longitude)

    # Convert heights to array
    heights_arr: np.ndarray
    if heights is None:
        heights_arr = np.zeros(n_sites)
    else:
        heights_arr = np.atleast_1d(heights)
        if heights_arr.size == 1:
            heights_arr = np.full(n_sites, heights_arr[0])
        elif heights_arr.size != n_sites:
            msg = (
                "Heights array must have same length as coordinates or be a single"
                " value"
            )
            raise ValueError(msg)

    # Convert rotation rates from deg/Myr to rad/year for internal calculations
    # deg/Myr × (π rad/180°) × (1 Myr/1e6 year) = rad/year  # noqa: RUF003
    W = np.array(
        [
            np.radians(wx) / MYR_TO_YEAR,  # rad/year
            np.radians(wy) / MYR_TO_YEAR,  # rad/year
            np.radians(wz) / MYR_TO_YEAR,  # rad/year
        ]
    )

    # Calculate modeled velocities for each site
    modeled_velocities = np.zeros((n_sites, 2))

    # Calculate uncertainties if covariance is provided
    if rotation_covariance is not None:
        uncertainties = np.zeros((n_sites, 2))

    for i, (lon, lat, hgt) in enumerate(
        zip(longitude, latitude, heights_arr, strict=False)
    ):
        # Get local frame matrix for this site
        local_frame = get_local_frame(lon, lat, hgt)

        # Calculate modeled velocity: v = R × ω  # noqa: RUF003
        # local_frame is 3x3, W is 3x1, result is 3x1 [ve, vn, vu]
        velocity_vector = local_frame @ W

        # Extract east and north components (first two elements)
        # Convert from m/year to mm/year
        modeled_velocities[i, 0] = velocity_vector[0] * 1000  # East component (mm/year)
        modeled_velocities[i, 1] = (
            velocity_vector[1] * 1000
        )  # North component (mm/year)

        # Calculate uncertainties if covariance is provided
        if rotation_covariance is not None:
            # Extract East-North components (2x3 matrix)
            J = local_frame[:2, :]  # Jacobian matrix for this site

            # Propagate uncertainty: Σv = J × Σω × J^T  # noqa: RUF003
            # where Σω is the rotation parameter covariance matrix
            velocity_covariance = J @ rotation_covariance @ J.T

            # Extract standard deviations and convert from m/year to mm/year
            uncertainties[i, 0] = (
                np.sqrt(velocity_covariance[0, 0]) * 1000
            )  # σve (mm/year)  # noqa: RUF003
            uncertainties[i, 1] = (
                np.sqrt(velocity_covariance[1, 1]) * 1000
            )  # σvn (mm/year)  # noqa: RUF003

    # Return results based on whether uncertainties were requested
    if rotation_covariance is None:
        # Return single pair for single input, array for multiple inputs
        if n_sites == 1:
            return modeled_velocities[0]
        else:
            return modeled_velocities
    else:
        # Always return arrays when uncertainties are included
        if n_sites == 1:
            return modeled_velocities[0], uncertainties[0]
        else:
            return modeled_velocities, uncertainties


def model_velocities_from_euler_pole(
    longitude: float | np.ndarray,
    latitude: float | np.ndarray,
    euler_longitude: float,
    euler_latitude: float,
    omega: float,
    heights: float | np.ndarray | None = None,
    euler_covariance: np.ndarray | None = None,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Compute plate motion velocities from Euler pole parameters.

    Parameters
    ----------
    longitude : float or ndarray
        Site longitude(s) in degrees.
    latitude : float or ndarray
        Site latitude(s) in degrees.
    euler_longitude : float
        Euler pole longitude in degrees.
    euler_latitude : float
        Euler pole latitude in degrees.
    omega : float
        Angular velocity in deg/Myr.
    heights : float or ndarray, optional
        Site heights above ellipsoid in meters. Single value applied to all
        sites if provided. Default is None (0 m).
    euler_covariance : ndarray, optional
        3x3 covariance matrix of Euler pole rotation parameters in (rad/year)^2.
        If provided, velocity uncertainties are returned.

    Returns
    -------
    ndarray or tuple of ndarray
        If `euler_covariance` is None:
            Array of modeled velocities [ve, vn] in mm/year
            (shape: (n_sites, 2) or (2,) for single site).
        If `euler_covariance` is provided:
            Tuple of (modeled_velocities, velocity_uncertainties)
            - modeled_velocities : ndarray of shape (n_sites, 2) in mm/year
            - velocity_uncertainties : ndarray of shape (n_sites, 2) in mm/year

    """
    # Convert Euler pole to rotation rate components
    wx, wy, wz = euler_pole_to_rotation_rate(euler_longitude, euler_latitude, omega)

    # Convert from rad/year back to deg/Myr for consistency with predict_plate_motion
    wx_deg_myr = np.degrees(wx) * MYR_TO_YEAR
    wy_deg_myr = np.degrees(wy) * MYR_TO_YEAR
    wz_deg_myr = np.degrees(wz) * MYR_TO_YEAR

    return model_plate_velocities(
        longitude,
        latitude,
        wx_deg_myr,
        wy_deg_myr,
        wz_deg_myr,
        heights,
        euler_covariance,
    )
