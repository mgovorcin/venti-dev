# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Sensor abstraction: everything sensor-specific in one place (PRD R-X1).

Algorithm code must not hard-code Sentinel-1 values. A `SensorSpec` carries the
radar wavelength (and hence the unwrapping cycle, half a wavelength of LOS
displacement), the product filename grammar, the dataset names of the layers
the calibration uses, and the static-layer grammar. DISP-S1 is implemented;
DISP-NISAR is registered as a known sensor whose spec lands with cal-disp v2
(plan T58), so asking for it fails with a clear message instead of silently
using C-band numbers.

Examples
--------
>>> spec = get_sensor("S1")
>>> round(spec.cycle_m, 5)
0.02773
>>> name = (
...     "OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z"
...     "_v1.0_20251027T005420Z.nc"
... )
>>> pid = spec.parse_filename(name)
>>> pid.frame_id, pid.reference_date.date().isoformat()
(8882, '2022-01-11')

"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import xarray as xr

logger = logging.getLogger(__name__)

__all__ = [
    "NISAR",
    "S1",
    "ProductId",
    "SensorSpec",
    "StaticLayerId",
    "UnsupportedSensorError",
    "get_sensor",
    "sensor_for_file",
]

_TIMESTAMP = "%Y%m%dT%H%M%SZ"


class UnsupportedSensorError(ValueError):
    """The file or sensor name is not one Venti knows how to read."""


@dataclass(frozen=True)
class ProductId:
    """What a displacement product's filename says about it."""

    path: Path
    sensor: str
    frame_id: int
    reference_date: datetime
    secondary_date: datetime
    polarization: str
    version: str
    production_date: datetime
    mode: str = "IW"

    @property
    def span_days(self) -> float:
        """Temporal baseline of the pair in days."""
        return (self.secondary_date - self.reference_date).total_seconds() / 86400.0


@dataclass(frozen=True)
class StaticLayerId:
    """What a static-layer filename says about it."""

    path: Path
    sensor: str
    frame_id: int
    reference_date: datetime
    platform: str
    version: str
    layer: str


@dataclass(frozen=True)
class SensorSpec:
    """Sensor-specific constants and readers for one DISP product family.

    Attributes
    ----------
    name : str
        Short name used in configs and metadata (``"S1"``, ``"NISAR"``).
    wavelength_m : float
        Nominal one-way radar wavelength in metres. The value written in a
        product's ``/identification/radar_wavelength`` takes precedence when a
        file is at hand (`read_wavelength`).
    product_pattern : re.Pattern
        Filename grammar of the displacement product.
    static_pattern : re.Pattern
        Filename grammar of the static layers (DEM, LOS, ...).
    displacement_layer, mask_layers, quality_layers : str / tuple[str, ...]
        Dataset names in the product's main group.
    corrections_group, correction_layers : str / tuple[str, ...]
        Group holding the correction layers and the names Venti knows.
    static_layers : tuple[str, ...]
        Static-layer names.
    implemented : bool
        False for sensors that are registered but whose readers are not
        written yet; `get_sensor` refuses them.

    """

    name: str
    wavelength_m: float
    product_pattern: re.Pattern[str]
    static_pattern: re.Pattern[str]
    displacement_layer: str = "displacement"
    mask_layers: tuple[str, ...] = ("recommended_mask", "water_mask")
    quality_layers: tuple[str, ...] = (
        "temporal_coherence",
        "timeseries_inversion_residuals",
        "connected_component_labels",
    )
    corrections_group: str = "corrections"
    correction_layers: tuple[str, ...] = (
        "solid_earth_tide",
        "ionospheric_delay",
        "perpendicular_baseline",
    )
    static_layers: tuple[str, ...] = ("dem", "line_of_sight_enu", "layover_shadow_mask")
    implemented: bool = True
    notes: str = field(default="", compare=False)

    # -- constants -------------------------------------------------------------

    @property
    def cycle_m(self) -> float:
        """One unwrapping cycle of two-way LOS displacement: half the wavelength."""
        return self.wavelength_m / 2.0

    # -- filenames -------------------------------------------------------------

    def matches(self, path: Path | str) -> bool:
        """Return True if `path` has this sensor's displacement-product filename."""
        return self.product_pattern.match(Path(path).name) is not None

    def parse_filename(self, path: Path | str) -> ProductId:
        """Parse a displacement product filename into a `ProductId`."""
        path = Path(path)
        match = self.product_pattern.match(path.name)
        if match is None:
            msg = f"{path.name} is not a {self.name} displacement product filename"
            raise UnsupportedSensorError(msg)
        g = match.groupdict()
        return ProductId(
            path=path,
            sensor=self.name,
            frame_id=int(g["frame_id"]),
            reference_date=datetime.strptime(g["reference"], _TIMESTAMP),
            secondary_date=datetime.strptime(g["secondary"], _TIMESTAMP),
            polarization=g["pol"],
            version=g["version"],
            production_date=datetime.strptime(g["production"], _TIMESTAMP),
            mode=g.get("mode", "IW"),
        )

    def parse_static_filename(self, path: Path | str) -> StaticLayerId:
        """Parse a static-layer filename into a `StaticLayerId`."""
        path = Path(path)
        match = self.static_pattern.match(path.name)
        if match is None:
            msg = f"{path.name} is not a {self.name} static-layer filename"
            raise UnsupportedSensorError(msg)
        g = match.groupdict()
        return StaticLayerId(
            path=path,
            sensor=self.name,
            frame_id=int(g["frame_id"]),
            reference_date=datetime.strptime(g["date"], "%Y%m%d"),
            platform=g["platform"],
            version=g["version"],
            layer=g["layer"],
        )

    # -- readers -----------------------------------------------------------------

    def read_wavelength(self, path: Path | str) -> float:
        """Radar wavelength in metres from the product, else the nominal value.

        Reads ``/identification/radar_wavelength``; a product without it
        (synthetic test files, older versions) falls back to `wavelength_m`
        with a warning, so the sensor default is never applied silently to a
        file that disagrees with it.
        """
        import h5py

        with h5py.File(path, "r") as f:
            node = f.get("identification/radar_wavelength")
            if node is not None:
                value = float(np.asarray(node[()]).item())
                if abs(value - self.wavelength_m) / self.wavelength_m > 0.05:
                    msg = (
                        f"{Path(path).name} says radar_wavelength={value:.5f} m, "
                        f"far from the {self.name} nominal {self.wavelength_m:.5f} m"
                    )
                    raise UnsupportedSensorError(msg)
                return value
        logger.warning(
            "%s has no /identification/radar_wavelength; using the %s nominal %.5f m",
            Path(path).name,
            self.name,
            self.wavelength_m,
        )
        return self.wavelength_m

    def open_displacement(self, path: Path | str) -> xr.DataArray:
        """Open the displacement layer (lazily) as an `xarray.DataArray`."""
        import xarray as xr

        ds = xr.open_dataset(path, engine="h5netcdf")
        return ds[self.displacement_layer]

    def open_mask(self, path: Path | str, layer: str) -> np.ndarray:
        """Read a mask layer as a boolean array (True = valid)."""
        import xarray as xr

        if layer not in self.mask_layers:
            msg = f"{layer!r} is not a {self.name} mask layer {self.mask_layers}"
            raise ValueError(msg)
        with xr.open_dataset(path, engine="h5netcdf") as ds:
            if layer not in ds:
                msg = f"{Path(path).name} has no {layer!r} layer"
                raise KeyError(msg)
            return ds[layer].to_numpy().astype(bool)

    def available_corrections(self, path: Path | str) -> tuple[str, ...]:
        """Names of the known correction layers present in the product."""
        import h5py

        with h5py.File(path, "r") as f:
            group = f.get(self.corrections_group)
            if group is None:
                return ()
            return tuple(name for name in self.correction_layers if name in group)

    def open_correction(self, path: Path | str, layer: str) -> np.ndarray | None:
        """Read a correction layer, or None if the product does not carry it."""
        from .io import read_netcdf_correction

        if layer not in self.correction_layers:
            msg = (
                f"{layer!r} is not a {self.name} correction layer "
                f"{self.correction_layers}"
            )
            raise ValueError(msg)
        return read_netcdf_correction(path, layer, group=self.corrections_group)


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------

S1 = SensorSpec(
    name="S1",
    wavelength_m=0.05546576,  # Sentinel-1 C-band, 5.405 GHz
    product_pattern=re.compile(
        r"OPERA_L3_DISP-S1_(?P<mode>\w+)_F(?P<frame_id>\d+)_(?P<pol>\w+)_"
        r"(?P<reference>\d{8}T\d{6}Z)_(?P<secondary>\d{8}T\d{6}Z)_"
        r"v(?P<version>[\d.]+)_(?P<production>\d{8}T\d{6}Z)\.nc$"
    ),
    static_pattern=re.compile(
        r"OPERA_L3_DISP-S1-STATIC_F(?P<frame_id>\d+)_(?P<date>\d{8})_"
        r"(?P<platform>S1[ABCD])_v(?P<version>[\d.]+)_(?P<layer>[\w]+)\.tif$"
    ),
)

NISAR = SensorSpec(
    name="NISAR",
    wavelength_m=0.2413,  # L-band, 1.2425 GHz; one cycle ~ 12 cm
    product_pattern=re.compile(r"OPERA_L3_DISP-NISAR_.*\.(nc|h5)$"),
    static_pattern=re.compile(r"OPERA_L3_DISP-NISAR-STATIC_.*\.tif$"),
    correction_layers=("solid_earth_tide", "ionosphere_phase_screen"),
    implemented=False,
    notes=(
        "DISP-NISAR readers land with cal-disp v2 (plan T58); the displacement "
        "arrives ionosphere-corrected and the phase screen is a correction layer."
    ),
)

_REGISTRY: dict[str, SensorSpec] = {S1.name: S1, NISAR.name: NISAR}


def get_sensor(name: str) -> SensorSpec:
    """Return the implemented `SensorSpec` called `name` (case-insensitive)."""
    spec = _REGISTRY.get(name.upper())
    if spec is None:
        msg = f"Unknown sensor {name!r}; known: {sorted(_REGISTRY)}"
        raise UnsupportedSensorError(msg)
    if not spec.implemented:
        msg = f"Sensor {spec.name} is registered but not implemented yet. {spec.notes}"
        raise UnsupportedSensorError(msg)
    return spec


def sensor_for_file(path: Path | str) -> SensorSpec:
    """Pick the implemented sensor whose product grammar matches `path`."""
    name = Path(path).name
    for spec in _REGISTRY.values():
        if spec.product_pattern.match(name):
            if not spec.implemented:
                msg = f"{name} is a {spec.name} product; {spec.notes}"
                raise UnsupportedSensorError(msg)
            return spec
    msg = f"{name} does not match any known displacement product filename"
    raise UnsupportedSensorError(msg)
