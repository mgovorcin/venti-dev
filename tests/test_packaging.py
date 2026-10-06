# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
"""Regression tests for the packaging metadata in pyproject.toml.

pixi >= 0.48 refuses a manifest that declares the same package both as a PEP 621
dependency and as a `[tool.pixi.*pypi-dependencies]` entry ("X is already a
dependency"). `dem-stitcher` and `geepers` were declared twice on 2026-10-05;
these tests keep that from coming back, and pin the shared task/environment
names the other repos and CI rely on.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest

PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


@pytest.fixture(scope="module")
def pyproject() -> dict:
    with PYPROJECT.open("rb") as f:
        return tomllib.load(f)


def _name(requirement: str) -> str:
    """Normalise a PEP 508 requirement to its distribution name."""
    name = re.split(r"[\s\[<>=!~;@]", requirement.strip(), maxsplit=1)[0]
    return name.lower().replace("_", "-")


def _pep621_names(pyproject: dict) -> set[str]:
    project = pyproject["project"]
    names = {_name(r) for r in project.get("dependencies", [])}
    for extra in project.get("optional-dependencies", {}).values():
        names |= {_name(r) for r in extra}
    return names


def _pixi_pypi_names(pyproject: dict) -> set[str]:
    pixi = pyproject["tool"]["pixi"]
    names = {_name(k) for k in pixi.get("pypi-dependencies", {})}
    for feature in pixi.get("feature", {}).values():
        names |= {_name(k) for k in feature.get("pypi-dependencies", {})}
    return names


def test_no_package_declared_twice_for_pixi(pyproject):
    duplicates = _pep621_names(pyproject) & _pixi_pypi_names(pyproject)
    # The project itself is installed editable through pixi; that is the one
    # allowed overlap.
    duplicates.discard("venti")
    assert (
        not duplicates
    ), f"declared both as PEP 621 and pixi pypi dependency: {sorted(duplicates)}"


def test_shared_pixi_task_names_exist(pyproject):
    tasks = pyproject["tool"]["pixi"]["tasks"]
    assert {"test", "lint", "docs", "golden", "e2e"} <= set(tasks)


def test_shared_pixi_environments_exist(pyproject):
    envs = pyproject["tool"]["pixi"]["environments"]
    assert "dev" in envs
    assert "ops" in envs
    # `ops` is the operational set: no extra features.
    assert envs["ops"]["features"] == []
