#!/usr/bin/env python
# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
r"""Benchmark calibration variants on the F08882 golden pair at real GNSS stations.

Plan T30.4 (and the station metric of the trade study
``surface_fitting/gnss_pair_validation.py``): calibrate the golden pair
(2022-01-11 -> 2022-07-22) with `venti.calibration.two_pass.calibrate_pair`
for several option sets, then compare ``DISP - calibration`` with 123 real
UNR stations (IGS20 daily positions, not the grid the surface is fitted
to). InSAR per station = median of the valid pixels in a 7 x 7 window;
GNSS reference = "pair" (mean position +/-3 d around each date, n = 83) and
"midas" (MIDAS rate x 0.526 yr, n = 123). The station table is
``tests/data/f08882_golden_pair_stations.csv`` (built by the trade study).

Inputs: the cal-disp golden directory (DISP product, LOS static layer, the
staged UNR grid nodes + lookup). The GNSS LOS field is sampled with
`venti.gnss.sampling` at the fit-grid centres (as cal-disp does) unless
``--gnss-los`` points at a full-resolution ``.npy`` (metres).

Example::

    NUMPY_MADVISE_HUGEPAGE=0 python scripts/benchmark_golden_pair.py \\
        --golden-dir ../cal-disp/test_golden --out docs/benchmarks_f08882.md

"""

from __future__ import annotations

import argparse
import json
import logging
import resource
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import rasterio
import yaml

from venti.calibration.two_pass import CalibrationResult, calibrate_pair
from venti.gnss.sampling import GnssGridConfig, project_field_to_los, sample_gnss_enu
from venti.spatial.resample import upsample_array
from venti.workflow.config import AlgorithmParameters, CalibrationOptions

logger = logging.getLogger("benchmark")

HERE = Path(__file__).resolve().parent
STATIONS = HERE.parent / "tests" / "data" / "f08882_golden_pair_stations.csv"
UTM_EPSG = 32615
WIN = 7
MIN_PIX = 10
DT_YEARS = (pd.Timestamp("2022-07-22") - pd.Timestamp("2022-01-11")).days / 365.25


def variants(
    gamma_params: AlgorithmParameters, with_defo: bool
) -> dict[str, tuple[CalibrationOptions, bool]]:
    """Return ``name -> (options, exclude defo areas)``; tropo off, SET on, factor 6."""
    gamma = gamma_params.calibration_options.model_copy(deep=True)
    gamma.apply_tropo_correction = False

    def loclin(**weights: Any) -> CalibrationOptions:
        o = CalibrationOptions(
            unwrap_error_correction=False,
            apply_tropo_correction=False,
            downsample_factor=6,
            surface={
                "method": "loclin",
                "cutoff_wavelength_meters": 50_000.0,
                "fill_gaps": True,
                "two_pass": True,
            },
            weights={"coherence_power": 0.0, "robust": False, **weights},
            uncertainty={"k_grid": 3.9},
        )
        return o

    out: dict[str, tuple[CalibrationOptions, bool]] = {
        "gamma 600 km (golden config, tropo off)": (gamma, False),
        "loclin 50 km, two-pass, unweighted": (loclin(), False),
        "loclin 50 km, two-pass, coh^8": (loclin(coherence_power=8.0), False),
        "loclin 50 km, two-pass, coh^8 + robust": (
            loclin(coherence_power=8.0, robust=True),
            False,
        ),
    }
    if with_defo:
        out["loclin 50 km, two-pass, unweighted, defo excluded"] = (loclin(), True)
        out["loclin 50 km, two-pass, coh^8 + robust, defo excluded"] = (
            loclin(coherence_power=8.0, robust=True),
            True,
        )
    return out


def defo_mask(path: Path, data: dict[str, Any]) -> np.ndarray:
    """Remove-restore mask from a curated defo-area GeoJSON (True = excluded).

    The trade-study files predate the ``version`` member `AreaDB` requires;
    a missing version is recorded as the file name.
    """
    from venti.calibration.remove_restore import AreaDB, remove_restore_mask

    with path.open() as fh:
        geo = json.load(fh)
    geo.setdefault("version", f"file:{path.name}")
    areas = AreaDB.from_geojson(geo)
    mask = remove_restore_mask((data["x"], data["y"]), UTM_EPSG, areas)
    logger.info(
        "defo areas %s (%s): %d pixels (%.1f %%) excluded from the fit",
        [a.id for a in areas.items],
        areas.version,
        int(mask.sum()),
        100.0 * mask.mean(),
    )
    return mask


def load_disp(path: Path) -> dict[str, Any]:
    """Read the layers the benchmark needs from the DISP product."""
    with h5py.File(path) as h:
        disp = h["displacement"][()].astype(np.float32)
        out = {
            "disp": disp,
            "recommended": h["recommended_mask"][()].astype(bool),
            "water": h["water_mask"][()].astype(bool),
            "coherence": h["temporal_coherence"][()].astype(np.float32),
            "set": h["corrections/solid_earth_tide"][()].astype(np.float32),
            "x": h["x"][()],
            "y": h["y"][()],
            "wavelength": float(h["identification/radar_wavelength"][()]),
        }
        out["crs_wkt"] = h["spatial_ref"].attrs.get("crs_wkt")
    out["ref_point"] = reference_point(out)
    return out


def reference_point(data: dict[str, Any]) -> tuple[int, int]:
    """Choose the reference pixel the way cal-disp does.

    ``opera_utils`` ``find_reference_point`` on the temporal coherence over
    valid pixels (the product's own reference pixel can be NaN in a pair);
    falls back to the best pixel.
    """
    import tempfile

    mask = data["recommended"] & data["water"] & np.isfinite(data["disp"])
    quality = np.where(mask & np.isfinite(data["coherence"]), data["coherence"], 0.0)
    quality = quality.astype(np.float32)
    try:
        from opera_utils.disp import rebase_reference
    except ImportError:
        rebase_reference = None
    if rebase_reference is not None:
        x, y = data["x"], data["y"]
        dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
        with tempfile.TemporaryDirectory() as tmp:
            qf = Path(tmp) / "reference_quality.tif"
            with rasterio.open(
                qf,
                "w",
                driver="GTiff",
                height=quality.shape[0],
                width=quality.shape[1],
                count=1,
                dtype="float32",
                crs=data["crs_wkt"],
                transform=rasterio.transform.from_origin(
                    float(x[0]) - dx / 2, float(y[0]) - dy / 2, dx, -dy
                ),
            ) as dst:
                dst.write(quality, 1)
            try:
                row, col = rebase_reference.find_reference_point(qf)
                return int(row), int(col)
            except ValueError:
                pass
    row, col = np.unravel_index(int(np.argmax(quality)), quality.shape)
    return int(row), int(col)


def gnss_los_field(
    golden: Path, x: np.ndarray, y: np.ndarray, factor: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return GNSS LOS displacement and sigma (m) on the full grid.

    Sampled at the fit-grid centres, as cal-disp does, then upsampled.
    """
    gnss_dir = golden / "input_data" / "gnss"
    lookups = sorted(gnss_dir.glob("grid_latlon_lookup*.txt"))
    if not lookups:
        msg = f"no grid lookup in {gnss_dir}"
        raise FileNotFoundError(msg)
    cfg = GnssGridConfig(
        grid_lookup=lookups[0],
        station_dir=gnss_dir,
        utm_epsg=UTM_EPSG,
        buffer_meters=0.0,
    )
    c = factor // 2
    xc, yc = x[c::factor], y[c::factor]
    field = sample_gnss_enu(cfg, (xc, yc))
    los_file = next(
        (golden / "input_data" / "static_input").glob("*line_of_sight_enu.tif")
    )
    with rasterio.open(los_file) as ds:
        los = ds.read().astype(np.float64)
    e, n, u = (band[c::factor, c::factor][: len(yc), : len(xc)] for band in los)
    los_mm, sigma_mm = project_field_to_los(field, e, n, u, DT_YEARS)
    logger.info(
        "GNSS field: %d nodes, %d missing, digest %s",
        field.provenance.n_nodes,
        field.provenance.n_missing,
        field.provenance.digest[:12],
    )
    shape = (len(y), len(x))
    full = upsample_array(np.asarray(los_mm, dtype=np.float64), shape)
    full_sigma = upsample_array(
        np.nan_to_num(np.asarray(sigma_mm, dtype=np.float64), nan=0.0), shape
    )
    return (full / 1000.0).astype(np.float32), (full_sigma / 1000.0).astype(np.float32)


def station_stats(
    calibrated: np.ndarray, valid: np.ndarray, stations: pd.DataFrame
) -> dict[str, dict[str, float]]:
    """Residual statistics of the calibrated DISP against the stations (mm, LOS)."""
    h = WIN // 2
    insar = np.full(len(stations), np.nan)
    for i, (r, c) in enumerate(zip(stations["row"], stations["col"], strict=True)):
        win = (slice(int(r) - h, int(r) + h + 1), slice(int(c) - h, int(c) + h + 1))
        m = valid[win]
        if m.sum() >= MIN_PIX:
            insar[i] = float(np.nanmedian(calibrated[win][m])) * 1000.0
    out = {}
    for ref in ("pair", "midas"):
        g = stations[f"gnss_los_{ref}_mm"].to_numpy()
        ok = np.isfinite(g) & np.isfinite(insar)
        d = insar[ok] - g[ok]
        out[ref] = {
            "n": int(ok.sum()),
            "bias": float(d.mean()),
            "median": float(np.median(d)),
            "rmse": float(np.sqrt(np.mean(d**2))),
            "std": float(d.std()),
            "nmad": float(1.4826 * np.median(np.abs(d - np.median(d)))),
            "corr": (
                float(np.corrcoef(insar[ok], g[ok])[0, 1])
                if ok.sum() > 2
                else float("nan")
            ),
        }
    return out


def grid_check(gnss_los: np.ndarray, stations: pd.DataFrame) -> dict[str, float]:
    """Return UNR grid LOS minus the stations at the station pixels (grid fidelity)."""
    g = gnss_los[stations["row"].to_numpy(int), stations["col"].to_numpy(int)] * 1000.0
    out = {}
    for ref in ("pair", "midas"):
        s = stations[f"gnss_los_{ref}_mm"].to_numpy()
        ok = np.isfinite(s) & np.isfinite(g)
        d = g[ok] - s[ok]
        out[f"{ref}_bias"] = float(d.mean())
        out[f"{ref}_std"] = float(d.std())
        out[f"{ref}_n"] = int(ok.sum())
    return out


def run_variant(
    name: str,
    opts: CalibrationOptions,
    data: dict[str, Any],
    gnss_los: np.ndarray,
    gnss_los_std: np.ndarray | None,
    exclude: np.ndarray | None = None,
) -> tuple[CalibrationResult, float]:
    """Calibrate the pair with one option set; return the result and the seconds."""
    t0 = time.perf_counter()
    mask = data["recommended"] & data["water"] & np.isfinite(data["disp"])
    res = calibrate_pair(
        data["disp"],
        gnss_los,
        mask,
        data["ref_point"],
        opts,
        30.0,
        data["wavelength"] / 2.0,
        gnss_los_std=gnss_los_std,
        coherence=data["coherence"],
        exclude_mask=exclude,
        set_correction=data["set"],
        n_jobs=-1,
    )
    dt = time.perf_counter() - t0
    logger.info("%s: %.0f s, method %s, passes %d", name, dt, res.method, res.passes)
    return res, dt


def main(argv: list[str] | None = None) -> int:
    """Run the benchmark; return the exit code."""
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--golden-dir", type=Path, required=True)
    p.add_argument(
        "--gnss-los", type=Path, help="full-resolution GNSS LOS field (.npy, m)"
    )
    p.add_argument(
        "--defo-area",
        type=Path,
        help="defo-area GeoJSON (EPSG:4326); adds variants with these areas excluded",
    )
    p.add_argument("--out", type=Path, help="markdown report to write")
    p.add_argument("--json", type=Path, help="raw numbers")
    p.add_argument("--only", nargs="*", help="substrings of the variants to run")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    golden = args.golden_dir
    disp_file = next((golden / "input_data" / "disp").glob("*.nc"))
    data = load_disp(disp_file)
    logger.info(
        "DISP %s: %s, ref point %s",
        disp_file.name,
        data["disp"].shape,
        data["ref_point"],
    )
    with (golden / "configs" / "algorithm_parameters.yaml").open() as fh:
        gamma_params = AlgorithmParameters.from_dict(yaml.safe_load(fh))
    if args.gnss_los:
        gnss_los = np.load(args.gnss_los).astype(np.float32)
        gnss_los_std = None
    else:
        gnss_los, gnss_los_std = gnss_los_field(golden, data["x"], data["y"], 6)
    stations = pd.read_csv(STATIONS, index_col=0)
    check = grid_check(gnss_los, stations)
    logger.info("grid vs stations: %s", json.dumps(check))

    valid = data["recommended"] & np.isfinite(data["disp"])
    results: dict[str, Any] = {"grid_check": check, "variants": {}}
    raw = station_stats(data["disp"], valid, stations)
    results["variants"]["raw DISP (relative to its reference point)"] = {
        "stats": raw,
        "seconds": 0.0,
    }
    exclude = defo_mask(args.defo_area, data) if args.defo_area else None
    for name, (opts, use_defo) in variants(gamma_params, exclude is not None).items():
        if args.only and not any(s in name for s in args.only):
            continue
        res, dt = run_variant(
            name, opts, data, gnss_los, gnss_los_std, exclude if use_defo else None
        )
        calibrated = data["disp"] - res.calibration
        stats = station_stats(calibrated, valid, stations)
        entry = {
            "stats": stats,
            "seconds": dt,
            "method": res.method,
            "passes": res.passes,
            "fit_residual_std_mm": (
                None if res.fit_residual_std is None else 1000.0 * res.fit_residual_std
            ),
            "sigma_cal_median_mm": (
                None
                if res.sigma_cal is None
                else float(1000.0 * np.nanmedian(res.sigma_cal))
            ),
            "options": opts.model_dump(mode="json"),
            "defo_excluded": use_defo,
            "n_excluded_pixels": res.n_excluded_pixels,
        }
        results["variants"][name] = entry
        logger.info(
            "%s: pair RMSE %.1f (n=%d), midas RMSE %.1f (n=%d)",
            name,
            stats["pair"]["rmse"],
            stats["pair"]["n"],
            stats["midas"]["rmse"],
            stats["midas"]["n"],
        )
        del res, calibrated
    results["peak_rss_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

    lines = [
        (
            "| Variant | pair bias | pair RMSE | pair std | pair NMAD | midas bias |"
            " midas RMSE | midas std | midas NMAD | fit std | s |"
        ),
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name, e in results["variants"].items():
        s_p, s_m = e["stats"]["pair"], e["stats"]["midas"]
        fit = (
            ""
            if e.get("fit_residual_std_mm") is None
            else f"{e['fit_residual_std_mm']:.1f}"
        )
        lines.append(
            f"| {name} | {s_p['bias']:+.1f} | {s_p['rmse']:.1f} | {s_p['std']:.1f} |"
            f" {s_p['nmad']:.1f} | {s_m['bias']:+.1f} | {s_m['rmse']:.1f} |"
            f" {s_m['std']:.1f} | {s_m['nmad']:.1f} | {fit} | {e['seconds']:.0f} |"
        )
    n_pair = results["variants"]["raw DISP (relative to its reference point)"]["stats"][
        "pair"
    ]["n"]
    n_midas = results["variants"]["raw DISP (relative to its reference point)"][
        "stats"
    ]["midas"]["n"]
    header = (
        f"LOS, mm. Residual = InSAR - GNSS at the stations (pair n = {n_pair}, midas n"
        f" = {n_midas}). Grid vs stations: midas {check['midas_bias']:+.1f} +/-"
        f" {check['midas_std']:.1f} mm, pair {check['pair_bias']:+.1f} +/-"
        f" {check['pair_std']:.1f} mm. Peak RSS {results['peak_rss_gb']:.1f} GB."
    )
    report = "\n".join([header, "", *lines, ""])
    print(report)
    if args.out:
        args.out.write_text(report)
    if args.json:
        args.json.write_text(json.dumps(results, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
