# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Input staging for a DISP frame: DISP-S1 products, DEM, LOS geometry, tropo.

These were ``scripts/staging/*_cli.py`` and were imported through a
``sys.path`` insert; they are a package now so the stage workflow works from
an installed Venti. The scripts remain as thin wrappers. The staging
dependencies (asf_search, dem-stitcher, opera-utils[disp,asf]) are the
``research`` extra, not part of the calibration core.
"""
