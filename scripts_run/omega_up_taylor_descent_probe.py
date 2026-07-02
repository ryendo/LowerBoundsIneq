#!/usr/bin/env python3
"""Exploratory descent probe for the Taylor-side Omega_up certificate.

This script is intentionally separate from the production Omega_up/Omega_mid
verification code.  It tests the paper's Taylor-side idea for

    J_k,xx = 0.5*y*lambda_1,xx + R_k,xx

near the equilateral triangle.  The key diagnostic is whether the signed
first-order y Taylor model

    lambda_1,xx(x,y) >= lambda_1,xx(p0)
        + (-lambda_1,xxy(p0)) * (y0-y)

is numerically credible on a proposed rectangle.  If this signed model is used
pointwise in 0.5*y*lambda_1,xx + R_k,xx, the Taylor region can be pushed much
lower than the absolute-value third-derivative bound.

All computations here use double precision analytic equilateral eigenfunctions
and sampled Rxx values, so this is a preliminary diagnostic, not a proof
certificate.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import omega_up_taylor_probe as base  # noqa: E402


@dataclass(frozen=True)
class Rect:
    name: str
    rx: float
    ry: float


class FiniteSpectralModel:
    """Finite analytic-basis model for lambda_xx on nearby triangles."""

    def __init__(self, mode_count: int, nq: int) -> None:
        self.mode_count = mode_count
        self.nq = nq
        xq, yq, weights = base.triangle_quadrature(nq)
        modes = base.enumerate_equilateral_modes("D", mode_count)
        _, gx, gy = base.orthonormalize_by_clusters(modes, xq, yq, weights)
        self.gx = gx[:, :mode_count]
        self.gy = gy[:, :mode_count]
        weighted_gx = weights[:, None] * self.gx
        weighted_gy = weights[:, None] * self.gy
        self.gxx = self.gx.T @ weighted_gx
        self.gyy = self.gy.T @ weighted_gy
        self.gxy = self.gx.T @ weighted_gy + self.gy.T @ weighted_gx

    @staticmethod
    def _p_coefficients(x: float, y: float) -> tuple[float, float, float]:
        s = x - 0.5
        return 1.0 + s * s / (y * y), -s * base.Y0 / (y * y), base.Y0 * base.Y0 / (y * y)

    @staticmethod
    def _px_coefficients(x: float, y: float) -> tuple[float, float, float]:
        s = x - 0.5
        return 2.0 * s / (y * y), -base.Y0 / (y * y), 0.0

    @staticmethod
    def _pxx_coefficients(_x: float, y: float) -> tuple[float, float, float]:
        return 2.0 / (y * y), 0.0, 0.0

    def _matrix(self, coeffs: tuple[float, float, float]) -> np.ndarray:
        return coeffs[0] * self.gxx + coeffs[1] * self.gxy + coeffs[2] * self.gyy

    def lambda_xx(self, x: float, y: float) -> tuple[float, float]:
        """Finite spectral value of lambda_1,xx and lambda_1 at (x,y)."""
        a = self._matrix(self._p_coefficients(x, y))
        ax = self._matrix(self._px_coefficients(x, y))
        axx = self._matrix(self._pxx_coefficients(x, y))
        a = 0.5 * (a + a.T)
        ax = 0.5 * (ax + ax.T)
        axx = 0.5 * (axx + axx.T)
        eigvals, eigvecs = np.linalg.eigh(a)
        u = eigvecs[:, 0]
        b = eigvecs.T @ (ax @ u)
        lxx = float(u @ (axx @ u) - 2.0 * np.sum((b[1:] ** 2) / (eigvals[1:] - eigvals[0])))
        return lxx, float(eigvals[0])


def rxx_lower_grid(kind: str, rx: float, ry: float, samples: int) -> float:
    return base.rxx_lower_sampled(kind, rx, ry, samples)


def coupled_linear_grid(
    rect: Rect,
    lambda_xx0: float,
    y_slope: float,
    samples: int,
) -> dict[str, object]:
    """Sample the signed y-linear Taylor lower model for Jxx on a rectangle."""
    xs = np.linspace(0.5, 0.5 + rect.rx, samples)
    ys = np.linspace(base.Y0 - rect.ry, base.Y0, samples)
    x_grid, y_grid = np.meshgrid(xs, ys, indexing="ij")
    u_grid = base.Y0 - y_grid
    lambda_lower = lambda_xx0 + y_slope * u_grid
    spectral = 0.5 * y_grid * lambda_lower
    j1 = spectral + base.rxx_value("J1", x_grid, y_grid)
    j2 = spectral + base.rxx_value("J2", x_grid, y_grid)
    i1 = np.unravel_index(np.argmin(j1), j1.shape)
    i2 = np.unravel_index(np.argmin(j2), j2.shape)
    return {
        "name": rect.name,
        "rx": rect.rx,
        "ry": rect.ry,
        "x_sup": 0.5 + rect.rx,
        "y_inf": base.Y0 - rect.ry,
        "lambda_xx0": lambda_xx0,
        "y_slope": y_slope,
        "min_J1": float(j1[i1]),
        "min_J1_x": float(x_grid[i1]),
        "min_J1_y": float(y_grid[i1]),
        "min_J1_lambda_lower": float(lambda_lower[i1]),
        "min_J2": float(j2[i2]),
        "min_J2_x": float(x_grid[i2]),
        "min_J2_y": float(y_grid[i2]),
        "min_J2_lambda_lower": float(lambda_lower[i2]),
        "samples": samples,
    }


def direct_residual_grid(
    model: FiniteSpectralModel,
    rect: Rect,
    lambda_xx0: float,
    y_slope: float,
    samples: int,
) -> dict[str, object]:
    """Sample finite lambda_xx against the signed y-linear Taylor model."""
    min_residual = math.inf
    min_residual_point: tuple[float, float, float, float] | None = None
    min_j1 = math.inf
    min_j1_point: tuple[float, float, float, float] | None = None
    for x in np.linspace(0.5, 0.5 + rect.rx, samples):
        for y in np.linspace(base.Y0 - rect.ry, base.Y0, samples):
            lxx, lam1 = model.lambda_xx(float(x), float(y))
            residual = lxx - (lambda_xx0 + y_slope * (base.Y0 - float(y)))
            j1 = 0.5 * float(y) * lxx + float(base.rxx_value("J1", np.array([[x]]), np.array([[y]]))[0, 0])
            if residual < min_residual:
                min_residual = residual
                min_residual_point = (float(x), float(y), lxx, lam1)
            if j1 < min_j1:
                min_j1 = j1
                min_j1_point = (float(x), float(y), lxx, lam1)
    assert min_residual_point is not None
    assert min_j1_point is not None
    return {
        "name": rect.name,
        "rx": rect.rx,
        "ry": rect.ry,
        "x_sup": 0.5 + rect.rx,
        "y_inf": base.Y0 - rect.ry,
        "min_residual": min_residual,
        "min_residual_x": min_residual_point[0],
        "min_residual_y": min_residual_point[1],
        "min_residual_lambda_xx": min_residual_point[2],
        "min_residual_lambda1": min_residual_point[3],
        "min_direct_J1": min_j1,
        "min_direct_J1_x": min_j1_point[0],
        "min_direct_J1_y": min_j1_point[1],
        "min_direct_J1_lambda_xx": min_j1_point[2],
        "min_direct_J1_lambda1": min_j1_point[3],
        "samples": samples,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode-count", type=int, default=220)
    parser.add_argument("--nq", type=int, default=180)
    parser.add_argument("--grid-samples", type=int, default=401)
    parser.add_argument("--direct-samples", type=int, default=31)
    parser.add_argument("--outdir", default="results/omega_up_taylor_descent_probe")
    args = parser.parse_args()

    outdir = REPO_ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    model = FiniteSpectralModel(args.mode_count, args.nq)
    lambda_xx0, _ = model.lambda_xx(0.5, base.Y0)
    # Finite-difference slope for lambda_xx(0.5,y0-u) at u=0.
    h = 5.0e-5
    lambda_down, _ = model.lambda_xx(0.5, base.Y0 - h)
    y_slope = (lambda_down - lambda_xx0) / h

    rects = [
        Rect("old_eps_0p122", 0.244, 0.122),
        Rect("coupled_linear_eps_0p152", 2.0 * 0.15206564769202943, 0.15206564769202943),
        Rect("mid_top_narrow", 0.0038733460138323716, base.Y0 - 0.744643343738),
        Rect("signed_prev_y_0p8472", 0.1878944815, base.Y0 - 0.8472),
    ]

    coupled_rows = [
        coupled_linear_grid(rect, lambda_xx0, y_slope, args.grid_samples)
        for rect in rects
    ]
    residual_rows = [
        direct_residual_grid(model, rect, lambda_xx0, y_slope, args.direct_samples)
        for rect in rects
    ]

    write_csv(outdir / "coupled_linear_taylor_grid.csv", coupled_rows)
    write_csv(outdir / "direct_residual_grid.csv", residual_rows)

    with (outdir / "summary.md").open("w") as f:
        f.write("# Omega_up Taylor Descent Probe\n\n")
        f.write("Exploratory double-precision diagnostic; not a proof certificate.\n\n")
        f.write(f"- finite basis modes: `{args.mode_count}`\n")
        f.write(f"- quadrature order: `{args.nq}`\n")
        f.write(f"- finite `lambda_xx(p0)`: `{lambda_xx0:.12g}`\n")
        f.write(f"- finite signed y slope `-lambda_xxy(p0)`: `{y_slope:.12g}`\n\n")
        f.write("## Coupled Signed Taylor Model\n\n")
        f.write("| region | rx | ry | y_inf | min J1 | min J2 | min J1 point |\n")
        f.write("|---|---:|---:|---:|---:|---:|---|\n")
        for row in coupled_rows:
            f.write(
                f"| {row['name']} | {row['rx']:.12g} | {row['ry']:.12g} | "
                f"{row['y_inf']:.12g} | {row['min_J1']:.8e} | {row['min_J2']:.8e} | "
                f"({row['min_J1_x']:.8g}, {row['min_J1_y']:.8g}) |\n"
            )
        f.write("\n## Direct Finite Residual Check\n\n")
        f.write("| region | min residual | residual point | min direct J1 | min direct J1 point |\n")
        f.write("|---|---:|---|---:|---|\n")
        for row in residual_rows:
            f.write(
                f"| {row['name']} | {row['min_residual']:.8e} | "
                f"({row['min_residual_x']:.8g}, {row['min_residual_y']:.8g}) | "
                f"{row['min_direct_J1']:.8e} | "
                f"({row['min_direct_J1_x']:.8g}, {row['min_direct_J1_y']:.8g}) |\n"
            )

    print(f"lambda_xx0={lambda_xx0:.12g}")
    print(f"y_slope={y_slope:.12g}")
    print(f"wrote {outdir}")


if __name__ == "__main__":
    main()
