# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""The staging CLIs are a package (plan T17): importable, and the scripts are shims."""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULES = [
    "dem_cli",
    "disp_cli",
    "los_cli",
    "tropo_cli",
    "stage_frame_cli",
    "stage_window_cli",
    "utils",
]


@pytest.mark.parametrize("name", MODULES)
def test_staging_module_imports(name):
    pytest.importorskip(
        "asf_search"
    )  # research extra; the CLIs import it at module level
    mod = importlib.import_module(f"venti.staging.{name}")
    assert mod.__name__ == f"venti.staging.{name}"


@pytest.mark.parametrize("name", [m for m in MODULES if m != "utils"])
def test_script_is_a_thin_wrapper(name):
    text = (ROOT / "scripts/staging" / f"{name}.py").read_text()
    assert f"from venti.staging.{name} import main" in text
    assert len(text.splitlines()) < 25


def test_stage_frame_data_has_no_sys_path_hack():
    text = (ROOT / "src/venti/workflow/stage_frame_data.py").read_text()
    assert "sys.path" not in text
    assert "_STAGING_DIR" not in text
    assert re.search(r"from venti\.staging\.dem_cli import", text)


def test_staging_utils_is_core_safe():
    """utils has no asf/dem-stitcher dependency, so the package itself imports in the core."""
    import venti.staging  # noqa: F401
    from venti.staging.utils import parse_date

    assert parse_date("2016-06-15").year == 2016
