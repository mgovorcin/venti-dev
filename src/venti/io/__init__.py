"""I/O utilities for venti.

This module provides I/O functionality for reading and writing
geospatial data formats:
- Functional API: read/write functions
- Object-oriented API: RasterReader and RasterWriter for high-level operations
"""

from __future__ import annotations

__all__ = [
    # Product classes
    "CalProduct",
    "NetCDFData",
    "RasterData",
    "RasterMetadata",
    # Object-oriented API
    "RasterReader",
    "RasterWriter",
    "VlmProduct",
    "get_bounds",
    # Functional API
    "read_geotiff",
    "read_netcdf",
    "read_netcdf_correction",
    "update_netcdf_variable",
    "write_geotiff",
]


def __getattr__(name: str):
    """Lazily import I/O functions and classes on first access."""
    # Object-oriented API - Read operations
    if name in ["RasterReader", "RasterData", "RasterMetadata", "NetCDFData"]:
        from .read import (
            NetCDFData,
            RasterData,
            RasterMetadata,
            RasterReader,
        )

        globals()["RasterReader"] = RasterReader
        globals()["RasterData"] = RasterData
        globals()["RasterMetadata"] = RasterMetadata
        globals()["NetCDFData"] = NetCDFData

        return globals()[name]

    # Object-oriented API - Write operations
    if name == "RasterWriter":
        from .write import RasterWriter

        globals()["RasterWriter"] = RasterWriter

        return RasterWriter

    # Product classes
    if name in ["CalProduct", "VlmProduct"]:
        from .product import CalProduct, VlmProduct

        globals()["CalProduct"] = CalProduct
        globals()["VlmProduct"] = VlmProduct

        return globals()[name]

    # Functional API
    if name in [
        "read_geotiff",
        "write_geotiff",
        "read_netcdf",
        "read_netcdf_correction",
        "update_netcdf_variable",
        "get_bounds",
    ]:
        from .raster import (
            get_bounds,
            read_geotiff,
            read_netcdf,
            read_netcdf_correction,
            update_netcdf_variable,
            write_geotiff,
        )

        # Cache imports
        globals()["read_geotiff"] = read_geotiff
        globals()["write_geotiff"] = write_geotiff
        globals()["read_netcdf"] = read_netcdf
        globals()["read_netcdf_correction"] = read_netcdf_correction
        globals()["update_netcdf_variable"] = update_netcdf_variable
        globals()["get_bounds"] = get_bounds

        return globals()[name]

    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
