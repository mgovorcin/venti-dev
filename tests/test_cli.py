# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Tests for the ``venti`` command-line interface, run as a real subprocess."""

from __future__ import annotations

import subprocess
import sys

import pytest


def _venti(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "venti", *args],
        capture_output=True,
        text=True,
        cwd=cwd,
        check=False,  # the tests inspect failing exit codes
    )


def test_config_writes_templates_that_load(tmp_path):
    from venti.workflow.config import AlgorithmParameters

    result = _venti("config", "--output-dir", str(tmp_path))

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "runconfig.yaml").exists()
    AlgorithmParameters.from_yaml(tmp_path / "algorithm_parameters.yaml")
    assert "Configuration templates created" in result.stderr


@pytest.mark.parametrize(
    "args",
    [
        ("run", "--config-file", "missing.yaml"),
        ("run-single", "--config-file", "missing.yaml", "--disp-file", "x.nc"),
    ],
)
def test_missing_config_fails_with_logged_error(tmp_path, args):
    result = _venti(*args, cwd=tmp_path)

    assert result.returncode == 1
    assert "ERROR" in result.stderr


def test_invalid_log_level_is_rejected(tmp_path):
    result = _venti(
        "run", "--config-file", "x.yaml", "--log-level", "LOUD", cwd=tmp_path
    )

    assert result.returncode != 0
    assert "LOUD" in result.stdout + result.stderr
