# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Frozen UNR grid snapshots (PRD R-G4; plan T48.1)."""

from __future__ import annotations

import importlib.util
import json
from datetime import date
from pathlib import Path

import numpy as np
import pytest

from venti.gnss.sampling import GnssGridConfig
from venti.gnss.snapshot import (
    grid_config_from_snapshot,
    load_snapshot,
    snapshot_unr_grid,
    verify_snapshot,
)
from venti.gnss.unr import read_station_rate

# a v0.3-style lookup: id, longitude in 0..360, latitude
LOOKUP = "\n".join(
    f"{i:06d} {lon:.6f} {lat:.6f}"
    for i, (lon, lat) in enumerate(
        [(264.5, 29.5), (264.7, 29.7), (265.0, 30.0), (10.0, 45.0), (200.0, -20.0)],
        start=1,
    )
)
FAILING = {3}


def fake_fetch(url: str) -> str:
    assert "grid_latlon_lookup.txt" in url
    return LOOKUP + "\n"


def fake_download(ids, out_dir: Path, grid_type: str) -> list[Path]:
    """Write geepers-style files ``<id>_IGS20.tenv8``; node 3 'fails'."""
    files = []
    for i in ids:
        if i in FAILING:
            continue
        years = np.arange(2015.0, 2025.0, 0.5)
        ve, vn, vu = -12.0 + i, 3.0, -1.5 * i  # mm/yr
        rows = np.c_[
            years,
            ve * (years - 2020),
            vn * (years - 2020),
            vu * (years - 2020),
            np.full(years.size, 0.3),
            np.full(years.size, 0.3),
            np.full(years.size, 0.9),
            np.zeros(years.size),
        ]
        f = out_dir / f"{i:06d}_IGS20.tenv8"
        np.savetxt(f, rows, fmt="%.6f")
        files.append(f)
    return files


@pytest.fixture
def snapshot(tmp_path):
    return snapshot_unr_grid(
        tmp_path,
        grid_types=("constant",),
        snapshot_date=date(2026, 10, 6),
        fetch_lookup=fake_fetch,
        download=fake_download,
        notes="test",
    )


def test_layout_manifest_and_info(snapshot):
    assert snapshot.name == "unr_grid_0.3_20261006"
    info = load_snapshot(snapshot)
    assert info.snapshot_id == "unr_grid_0.3_IGS20_20261006"
    assert info.unr_version == "0.3"
    assert info.reference_frame == "IGS20"
    assert info.grid_types == ["constant"]
    assert info.n_lookup_nodes == 5
    assert info.n_selected_nodes == 5
    assert info.n_downloaded == {"constant": 4}
    assert info.n_failed == {"constant": 1}
    assert info.failed_ids == {"constant": [3]}
    assert info.data_span["constant"] == [2015.0, 2024.5]
    assert info.bounds_snwe is None
    assert len(info.lookup_sha256) == 64
    assert "geodesy.unr.edu" in info.source["lookup_url"]
    assert info.notes == "test"
    # the lookup is stored byte for byte (0..360 longitudes untouched)
    assert (snapshot / "grid_latlon_lookup.txt").read_text() == LOOKUP + "\n"
    nodes = sorted(p.name for p in (snapshot / "nodes").iterdir())
    assert nodes == [f"{i:06d}_IGS20_constant.tenv8" for i in (1, 2, 4, 5)]
    manifest = (snapshot / "MANIFEST.sha256").read_text().splitlines()
    assert len(manifest) == 6  # lookup + snapshot.json + 4 nodes
    assert all(len(line.split()[0]) == 64 for line in manifest)
    assert verify_snapshot(snapshot) == []
    # no temp directory left behind
    assert not [p for p in snapshot.iterdir() if p.name.startswith(".dl_")]


def test_verify_detects_changes(snapshot):
    f = snapshot / "nodes" / "000001_IGS20_constant.tenv8"
    f.write_text(f.read_text() + "2025.0 0 0 0 0 0 0 0\n")
    (snapshot / "nodes" / "000002_IGS20_constant.tenv8").unlink()
    (snapshot / "extra.txt").write_text("x")
    problems = verify_snapshot(snapshot)
    assert sorted(problems) == [
        "changed nodes/000001_IGS20_constant.tenv8",
        "missing nodes/000002_IGS20_constant.tenv8",
        "unlisted extra.txt",
    ]
    assert verify_snapshot(snapshot.parent) == ["missing MANIFEST.sha256"]
    with pytest.raises(ValueError, match="fails verification"):
        grid_config_from_snapshot(snapshot, 32615, verify=True)


def test_bounds_select_nodes_and_both_grid_types(tmp_path):
    out = snapshot_unr_grid(
        tmp_path,
        grid_types=("constant", "variable"),
        bounds_snwe=(29.0, 31.0, -96.0, -94.0),  # Houston; lookup is in 0..360
        snapshot_date=date(2026, 1, 1),
        fetch_lookup=fake_fetch,
        download=fake_download,
    )
    info = load_snapshot(out)
    assert info.n_selected_nodes == 3
    assert info.bounds_snwe == [29.0, 31.0, -96.0, -94.0]
    assert info.n_downloaded == {"constant": 2, "variable": 2}
    assert info.failed_ids == {"constant": [3], "variable": [3]}
    names = sorted(p.name for p in (out / "nodes").iterdir())
    assert names == [
        "000001_IGS20_constant.tenv8",
        "000001_IGS20_variable.tenv8",
        "000002_IGS20_constant.tenv8",
        "000002_IGS20_variable.tenv8",
    ]
    with pytest.raises(FileExistsError, match="never overwritten"):
        snapshot_unr_grid(
            tmp_path,
            snapshot_date=date(2026, 1, 1),
            fetch_lookup=fake_fetch,
            download=fake_download,
        )
    with pytest.raises(ValueError, match="bounds"):
        snapshot_unr_grid(
            tmp_path,
            bounds_snwe=(31, 29, -96, -94),
            fetch_lookup=fake_fetch,
            download=fake_download,
        )
    with pytest.raises(ValueError, match="grid type"):
        snapshot_unr_grid(
            tmp_path,
            grid_types=("daily",),
            fetch_lookup=fake_fetch,
            download=fake_download,
        )

    # a failed build leaves nothing behind that could pass for a snapshot
    def boom(ids, out_dir, grid_type):
        msg = "network down"
        raise RuntimeError(msg)

    with pytest.raises(RuntimeError, match="network down"):
        snapshot_unr_grid(
            tmp_path,
            snapshot_date=date(2026, 2, 2),
            fetch_lookup=fake_fetch,
            download=boom,
        )
    assert not (tmp_path / "unr_grid_0.3_20260202").exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["unr_grid_0.3_20260101"]


def test_grid_config_from_snapshot_feeds_the_sampler(snapshot):
    cfg = grid_config_from_snapshot(
        snapshot, 32615, buffer_meters=50_000.0, verify=True
    )
    assert isinstance(cfg, GnssGridConfig)
    assert cfg.snapshot_id == "unr_grid_0.3_IGS20_20261006"
    assert cfg.version == "0.3"
    assert cfg.reference_frame == "IGS20"
    assert cfg.grid_type == "constant"
    assert cfg.buffer_meters == 50_000.0
    assert cfg.grid_lookup == snapshot / "grid_latlon_lookup.txt"
    assert cfg.station_file(1).exists()
    ve, vn, vu, se, sn, su = read_station_rate(cfg.station_file(1))
    assert (ve, vn, vu) == pytest.approx((-11.0, 3.0, -1.5), abs=1e-6)
    assert (se, sn, su) == (0.3, 0.3, 0.9)
    with pytest.raises(ValueError, match="not 'variable'"):
        grid_config_from_snapshot(snapshot, 32615, grid_type="variable")


def test_script_verify_and_argument_errors(snapshot, capsys):
    spec = importlib.util.spec_from_file_location(
        "snapshot_unr_grid",
        Path(__file__).parents[1] / "scripts" / "snapshot_unr_grid.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["--verify", str(snapshot)]) == 0
    assert "OK (4 node files)" in capsys.readouterr().out
    (snapshot / "nodes" / "000001_IGS20_constant.tenv8").write_text("broken\n")
    assert mod.main(["--verify", str(snapshot)]) == 1
    assert "changed nodes/000001" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        mod.main([])  # output_root required
    args = mod.build_parser().parse_args(
        [
            "/x",
            "--bounds",
            "24",
            "50",
            "-125",
            "-66",
            "--grid-types",
            "constant",
            "variable",
        ]
    )
    assert args.bounds == [24.0, 50.0, -125.0, -66.0]
    assert args.grid_types == [
        "constant",
        "variable",
    ]
    (snapshot / "snapshot.json").write_text(json.dumps({"bad": 1}))
    with pytest.raises(TypeError):
        load_snapshot(snapshot)
