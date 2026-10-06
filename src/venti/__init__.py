# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Venti package for calibrating OPERA DISP with GNSS."""

from typing import Any

from .log_setup import configure_logging

# Lazy imports to avoid pyproj initialization errors
_LAZY_MODULES = {
    "models": ".models",
    "unwrap": ".unwrap",
    "workflow": ".workflow",
    "gnss": ".gnss",
    "io": ".io",
    "spatial": ".spatial",
    "surface": ".surface",
}
# Main entry points available directly as ``venti.<name>``.
_LAZY_ATTRS = {
    "estimate_calibration_surface": ".surface",
    "CalibrationSurface": ".surface",
}


def __getattr__(name: str) -> Any:
    """Lazily import submodules on first access."""
    if name in _LAZY_MODULES:
        import importlib

        module = importlib.import_module(_LAZY_MODULES[name], package=__name__)
        # Cache the module
        globals()[name] = module
        return module

    if name in _LAZY_ATTRS:
        import importlib

        value = getattr(importlib.import_module(_LAZY_ATTRS[name], __name__), name)
        globals()[name] = value
        return value

    if name == "__version__":
        try:
            from ._version import __version__

            globals()["__version__"] = __version__
        except ImportError:
            return "unknown"
        else:
            return __version__

    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)


def __dir__():
    """List available attributes."""
    return [*_LAZY_MODULES, *_LAZY_ATTRS, "__version__", "configure_logging"]


__all__ = [
    "CalibrationSurface",
    "__version__",
    "configure_logging",
    "estimate_calibration_surface",
    "gnss",
    "io",
    "models",
    "spatial",
    "surface",
    "unwrap",
    "workflow",
]
