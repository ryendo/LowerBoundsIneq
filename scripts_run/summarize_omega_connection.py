#!/usr/bin/env python3
"""Summarize the Omega_mid / Omega_up connection.

This script does not run a new proof by itself.  It collects the certified
Omega_mid CSV output and the Omega_up Taylor/FEM diagnostic output and writes a
single report that says exactly which interface has been checked.

Paper formula being audited:

    J_{k,xx}(p) = 0.5 * y * lambda_{1,xx}(p) + R_{xx}^{(k)}(p).

The Taylor-side rows are exploratory unless their source explicitly says
otherwise; this is intentional, so that numerical evidence is not silently
promoted to a proof certificate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Y0 = math.sqrt(3.0) / 2.0
P0 = (0.5, Y0)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def float_or_nan(value: str | None) -> float:
    if value is None or value == "":
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def cell_geometry_stats(cell_file: Path) -> dict[str, Any]:
    rows = read_csv(cell_file)
    min_dist = math.inf
    closest: dict[str, Any] | None = None
    min_y = math.inf
    max_y = -math.inf
    min_x = math.inf
    max_x = -math.inf

    for row in rows:
        xs = [float_or_nan(row.get("x_inf")), float_or_nan(row.get("x_sup"))]
        ts = [float_or_nan(row.get("theta_inf")), float_or_nan(row.get("theta_sup"))]
        for x in xs:
            for theta in ts:
                if math.isnan(x) or math.isnan(theta):
                    continue
                y = x * math.tan(theta)
                dist = math.hypot(x - P0[0], y - P0[1])
                min_x = min(min_x, x)
                max_x = max(max_x, x)
                min_y = min(min_y, y)
                max_y = max(max_y, y)
                if dist < min_dist:
                    min_dist = dist
                    closest = {
                        "cell_id": row.get("i", ""),
                        "x": x,
                        "theta": theta,
                        "y": y,
                        "distance_to_p0": dist,
                    }

    return {
        "cell_file": str(cell_file),
        "cell_count": len(rows),
        "min_distance_to_p0": min_dist,
        "min_x": min_x,
        "max_x": max_x,
        "min_y": min_y,
        "max_y": max_y,
        "closest_corner": closest,
    }


def verified_result_stats(path: Path) -> dict[str, Any]:
    rows = read_csv(path)
    lower_values: list[float] = []
    verified = 0
    latest_by_cell: dict[str, dict[str, str]] = {}
    for row in rows:
        cid = row.get("cell_id")
        if cid:
            latest_by_cell[cid] = row
        flag = row.get("verified", row.get("is_verified", ""))
        if flag not in {"", "0", "false", "False"}:
            verified += 1
        val = float_or_nan(row.get("J_lower"))
        if not math.isnan(val):
            lower_values.append(val)

    latest_verified = 0
    latest_lowers: list[float] = []
    for row in latest_by_cell.values():
        flag = row.get("verified", "")
        if flag not in {"", "0", "false", "False"}:
            latest_verified += 1
        val = float_or_nan(row.get("J_lower"))
        if not math.isnan(val):
            latest_lowers.append(val)

    values = latest_lowers if latest_by_cell else lower_values
    return {
        "path": str(path),
        "row_count": len(rows),
        "unique_cell_count": len(latest_by_cell) if latest_by_cell else len(rows),
        "verified_rows": verified,
        "verified_latest_cells": latest_verified if latest_by_cell else verified,
        "min_J_lower": min(values) if values else math.nan,
        "max_J_lower": max(values) if values else math.nan,
    }


def omega_up_step_stats(path: Path) -> dict[str, Any]:
    rows = read_csv(path)
    vals: list[float] = []
    ok = 0
    inf_count = 0
    for row in rows:
        val = float_or_nan(row.get("L_lower"))
        if math.isinf(val):
            inf_count += 1
        elif not math.isnan(val):
            vals.append(val)
        if row.get("ok", "") not in {"", "0", "false", "False"}:
            ok += 1
    return {
        "path": str(path),
        "row_count": len(rows),
        "ok_rows": ok,
        "finite_min_L_lower": min(vals) if vals else math.nan,
        "finite_max_L_lower": max(vals) if vals else math.nan,
        "inf_rows": inf_count,
    }


def taylor_stats(taylor_dir: Path) -> dict[str, Any]:
    rows = read_csv(taylor_dir / "coupled_linear_taylor_grid.csv")
    residual_rows = read_csv(taylor_dir / "direct_residual_grid.csv")
    by_name = {row.get("name", ""): row for row in rows}
    residual_by_name = {row.get("name", ""): row for row in residual_rows}
    old = by_name.get("old_eps_0p122", {})
    old_res = residual_by_name.get("old_eps_0p122", {})
    return {
        "directory": str(taylor_dir),
        "rigor_status": "exploratory_double",
        "old_eps_0p122": {
            "rx": float_or_nan(old.get("rx")),
            "ry": float_or_nan(old.get("ry")),
            "min_J1": float_or_nan(old.get("min_J1")),
            "min_J2": float_or_nan(old.get("min_J2")),
            "lambda_xx0": float_or_nan(old.get("lambda_xx0")),
            "signed_y_slope": float_or_nan(old.get("y_slope")),
            "min_direct_residual": float_or_nan(old_res.get("min_residual")),
            "min_direct_J1": float_or_nan(old_res.get("min_direct_J1")),
        },
    }


def write_report(out_path: Path, data: dict[str, Any]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mid = data["omega_mid_geometry"]
    j1_mid = data["omega_mid_results"]["J1"]
    j2_mid = data["omega_mid_results"]["J2"]
    up_j1_x = data["omega_up_fem"]["J1_step_1_2"]
    up_j2_x = data["omega_up_fem"]["J2_step_1_2"]
    up_j1_y = data["omega_up_fem"]["J1_step_1_3"]
    up_j2_y = data["omega_up_fem"]["J2_step_1_3"]
    tay = data["omega_up_taylor"]["old_eps_0p122"]

    with out_path.open("w") as f:
        f.write("# Omega_mid / Omega_up Connection Summary\n\n")
        f.write("Generated by `scripts_run/summarize_omega_connection.py`.\n\n")
        f.write("## Interface Geometry\n\n")
        f.write(f"- Omega_mid cell file: `{mid['cell_file']}`\n")
        f.write(f"- Omega_mid cells: `{mid['cell_count']}`\n")
        f.write(f"- Closest Omega_mid cell corner to `p0`: `{mid['min_distance_to_p0']:.12g}`\n")
        f.write(f"- Claimed Omega_up epsilon: `{data['eps_up']:.12g}`\n")
        f.write(f"- Geometric overlap margin `eps_up - min_dist`: `{data['connection_margin']:.12g}`\n\n")

        f.write("## Omega_mid Certified CSVs\n\n")
        f.write("| functional | rows | unique cells | verified latest | min J_lower |\n")
        f.write("|---|---:|---:|---:|---:|\n")
        f.write(
            f"| J1 | {j1_mid['row_count']} | {j1_mid['unique_cell_count']} | "
            f"{j1_mid['verified_latest_cells']} | {j1_mid['min_J_lower']:.17e} |\n"
        )
        f.write(
            f"| J2 | {j2_mid['row_count']} | {j2_mid['unique_cell_count']} | "
            f"{j2_mid['verified_latest_cells']} | {j2_mid['min_J_lower']:.17e} |\n\n"
        )

        f.write("## Omega_up FEM/INTLAB Recalculation\n\n")
        f.write("| functional | step | rows | ok rows | finite min lower |\n")
        f.write("|---|---|---:|---:|---:|\n")
        f.write(f"| J1 | 1-2 x-convexity | {up_j1_x['row_count']} | {up_j1_x['ok_rows']} | {up_j1_x['finite_min_L_lower']:.17e} |\n")
        f.write(f"| J2 | 1-2 x-convexity | {up_j2_x['row_count']} | {up_j2_x['ok_rows']} | {up_j2_x['finite_min_L_lower']:.17e} |\n")
        f.write(f"| J1 | 1-3 y-axis | {up_j1_y['row_count']} | {up_j1_y['ok_rows']} | {up_j1_y['finite_min_L_lower']:.17e} |\n")
        f.write(f"| J2 | 1-3 y-axis | {up_j2_y['row_count']} | {up_j2_y['ok_rows']} | {up_j2_y['finite_min_L_lower']:.17e} |\n\n")

        f.write("## Omega_up Signed Taylor Diagnostic\n\n")
        f.write("This row evaluates the paper formula `J_{k,xx}=0.5*y*lambda_{xx}+R_{xx}^{(k)}` using the explicit equilateral spectral model and a signed y-linear Taylor lower model.  Its rigor flag is `exploratory_double`.\n\n")
        f.write("| region | rx | ry | lambda_xx(p0) | -lambda_xxy(p0) | min J1 | min J2 | direct residual min |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
        f.write(
            f"| old_eps_0p122 | {tay['rx']:.12g} | {tay['ry']:.12g} | "
            f"{tay['lambda_xx0']:.12g} | {tay['signed_y_slope']:.12g} | "
            f"{tay['min_J1']:.12g} | {tay['min_J2']:.12g} | {tay['min_direct_residual']:.12g} |\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--eps-up", type=float, default=0.122)
    parser.add_argument("--cell-file", default="inputs/cell_def.csv")
    parser.add_argument("--taylor-dir", default="results/omega_up_taylor_descent_probe")
    parser.add_argument("--out-md", default="results/omega_connection_summary.md")
    parser.add_argument("--out-json", default="results/omega_connection_summary.json")
    args = parser.parse_args()

    repo = args.repo.resolve()
    mid_geom = cell_geometry_stats(repo / args.cell_file)
    mid_geom["cell_file"] = args.cell_file
    data: dict[str, Any] = {
        "eps_up": args.eps_up,
        "omega_mid_geometry": mid_geom,
        "connection_margin": args.eps_up - mid_geom["min_distance_to_p0"],
        "omega_mid_results": {
            "J1": verified_result_stats(repo / "results" / "J1_OmegaMid.csv"),
            "J2": verified_result_stats(repo / "results" / "J2_OmegaMid.csv"),
        },
        "omega_up_fem": {
            "J1_step_1_2": omega_up_step_stats(repo / "results" / "J1_OmegaUp_step1_2_cells.csv"),
            "J2_step_1_2": omega_up_step_stats(repo / "results" / "J2_OmegaUp_step1_2_cells.csv"),
            "J1_step_1_3": omega_up_step_stats(repo / "results" / "J1_OmegaUp_step1_3_axis.csv"),
            "J2_step_1_3": omega_up_step_stats(repo / "results" / "J2_OmegaUp_step1_3_axis.csv"),
        },
        "omega_up_taylor": taylor_stats(repo / args.taylor_dir),
    }

    out_json = repo / args.out_json
    out_md = repo / args.out_md
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    write_report(out_md, data)
    print(f"Wrote {out_md}")
    print(f"Wrote {out_json}")


if __name__ == "__main__":
    main()
