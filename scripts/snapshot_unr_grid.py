#!/usr/bin/env python
# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
r"""Write or verify a frozen UNR grid snapshot (PRD R-G4; plan T48.1).

Examples
--------
::

    # the whole IGS20 constant grid, today's date
    python scripts/snapshot_unr_grid.py /data/unr_snapshots

    # constant + variable, CONUS only, for a reprocessing campaign
    python scripts/snapshot_unr_grid.py /data/unr_snapshots \\
        --grid-types constant variable --bounds 24 50 -125 -66 \\
        --notes "reprocessing 2026Q4"

    # check a snapshot against its manifest
    python scripts/snapshot_unr_grid.py --verify \\
        /data/unr_snapshots/unr_grid_0.3_20261006

"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import date
from pathlib import Path

from venti.gnss.snapshot import load_snapshot, snapshot_unr_grid, verify_snapshot


def build_parser() -> argparse.ArgumentParser:
    """Return the command-line parser."""
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "output_root", nargs="?", type=Path, help="parent directory for the snapshot"
    )
    p.add_argument("--version", default="0.3", help="UNR grid version (default 0.3)")
    p.add_argument("--frame", default="IGS20", help="reference frame (default IGS20)")
    p.add_argument(
        "--grid-types",
        nargs="+",
        default=["constant"],
        choices=["constant", "variable"],
        help="grid products to snapshot (default: constant)",
    )
    p.add_argument(
        "--bounds",
        nargs=4,
        type=float,
        metavar=("S", "N", "W", "E"),
        help="only nodes inside these degrees (default: the whole grid)",
    )
    p.add_argument(
        "--date", type=date.fromisoformat, help="snapshot date (default today)"
    )
    p.add_argument("--max-workers", type=int, default=8)
    p.add_argument(
        "--allow-missing",
        action="store_true",
        help="write the snapshot even if some nodes fail (ids recorded)",
    )
    p.add_argument("--notes", default="")
    p.add_argument(
        "--verify", type=Path, metavar="SNAPSHOT_DIR", help="verify instead of writing"
    )
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    """Write or verify a snapshot; return the exit code."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if args.verify is not None:
        problems = verify_snapshot(args.verify)
        info = load_snapshot(args.verify)
        if problems:
            print(f"{info.snapshot_id}: {len(problems)} problem(s)")
            for p in problems[:50]:
                print("  " + p)
            return 1
        print(f"{info.snapshot_id}: OK ({sum(info.n_downloaded.values())} node files)")
        return 0
    if args.output_root is None:
        build_parser().error("output_root is required unless --verify is given")
    out = snapshot_unr_grid(
        args.output_root,
        version=args.version,
        reference_frame=args.frame,
        grid_types=args.grid_types,
        bounds_snwe=args.bounds,
        snapshot_date=args.date,
        max_workers=args.max_workers,
        notes=args.notes,
        allow_missing=args.allow_missing,
    )
    info = load_snapshot(out)
    print(f"{info.snapshot_id} -> {out}")
    for gt in info.grid_types:
        print(f"  {gt}: {info.n_downloaded[gt]} nodes, {info.n_failed[gt]} failed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
