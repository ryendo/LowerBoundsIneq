#!/usr/bin/env python3
"""Exploratory Taylor spectral certificate near the equilateral triangle.

This script is intentionally isolated from the production Omega_up/Omega_mid
verification pipeline.  It estimates how large a near-equilateral rectangle can
be before the Taylor lower bound for J_xx loses positivity.

The computation uses analytic Lame-type equilateral triangle eigenfunctions,
but the spectral inner products and Rxx cell lower bounds are evaluated in
double precision.  Therefore every result produced here is marked
``exploratory_double`` and is not a proof certificate.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


Y0 = math.sqrt(3.0) / 2.0
AREA = math.sqrt(3.0) / 4.0
EIGEN_SCALE = 16.0 * math.pi**2 / 9.0


@dataclass(frozen=True)
class Mode:
    bc: str
    m: int
    n: int
    kind: str
    q: int
    eigenvalue: float
    label: str


@dataclass
class SpectralProbe:
    J: int
    lambda1: float
    lambda_x: float
    lambda_y: float
    lambda_xx0_lower: float
    m_xxx: float
    m_xxy: float
    tail_mode: str
    rigorous_flag: str


def parse_j_list(raw: str) -> list[int]:
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


def enumerate_equilateral_modes(bc: str, count: int) -> list[Mode]:
    """Return Lame-type analytic modes sorted by eigenvalue.

    Dirichlet modes use integer pairs 1 <= m <= n.  Positive Neumann modes use
    0 <= m <= n, excluding (0,0).  For m < n both symmetric and antisymmetric
    modes are included; for m = n the antisymmetric formula is identically zero.
    """
    limit = 1
    while True:
        modes: list[Mode] = []
        if bc == "D":
            for m in range(1, limit + 1):
                for n in range(m, limit + 1):
                    q = m * m + m * n + n * n
                    modes.append(Mode(bc, m, n, "s", q, EIGEN_SCALE * q, f"D({m},{n})s"))
                    if n > m:
                        modes.append(Mode(bc, m, n, "a", q, EIGEN_SCALE * q, f"D({m},{n})a"))
        elif bc == "N":
            for m in range(0, limit + 1):
                for n in range(m, limit + 1):
                    if m == 0 and n == 0:
                        continue
                    q = m * m + m * n + n * n
                    modes.append(Mode(bc, m, n, "s", q, EIGEN_SCALE * q, f"N({m},{n})s"))
                    if n > m:
                        modes.append(Mode(bc, m, n, "a", q, EIGEN_SCALE * q, f"N({m},{n})a"))
        else:
            raise ValueError(f"unknown boundary condition: {bc}")

        modes.sort(key=lambda z: (z.q, z.m, z.n, z.kind))
        if len(modes) >= count:
            return modes[:count]
        limit *= 2


def triangle_quadrature(nq: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tensor Gauss rule on the unit-side equilateral triangle."""
    z, w = np.polynomial.legendre.leggauss(nq)
    t = (z + 1.0) / 2.0
    wt = w / 2.0
    a = t[:, None]
    b = t[None, :]
    weights = (wt[:, None] * wt[None, :]) * Y0 * (1.0 - a)
    x = a + (1.0 - a) * b * 0.5
    y = (1.0 - a) * b * Y0
    return x.ravel(), y.ravel(), weights.ravel()


def raw_lame_mode(mode: Mode, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a raw Lame mode and its Cartesian gradient."""
    ell = -(mode.m + mode.n)
    u_coord = Y0 - y
    v_coord = math.sqrt(3.0) * (x - 0.5)
    triples = [
        (ell, mode.m - mode.n),
        (mode.m, mode.n - ell),
        (mode.n, ell - mode.m),
    ]

    value = np.zeros_like(x)
    ux = np.zeros_like(x)
    uy = np.zeros_like(x)
    for j, d in triples:
        a = math.pi * j / Y0
        b = math.pi * d / (3.0 * Y0)
        au = a * u_coord
        bv = b * v_coord
        if mode.bc == "D":
            sy = np.sin(au)
            cy = np.cos(au)
            if mode.kind == "s":
                cv = np.cos(bv)
                sv = np.sin(bv)
                value += sy * cv
                ux += sy * (-sv) * b * math.sqrt(3.0)
                uy += -a * cy * cv
            else:
                sv = np.sin(bv)
                cv = np.cos(bv)
                value += sy * sv
                ux += sy * cv * b * math.sqrt(3.0)
                uy += -a * cy * sv
        else:
            cy = np.cos(au)
            sy = np.sin(au)
            if mode.kind == "s":
                cv = np.cos(bv)
                sv = np.sin(bv)
                value += cy * cv
                ux += cy * (-sv) * b * math.sqrt(3.0)
                uy += a * sy * cv
            else:
                sv = np.sin(bv)
                cv = np.cos(bv)
                value += cy * sv
                ux += cy * cv * b * math.sqrt(3.0)
                uy += a * sy * sv
    return value, ux, uy


def orthonormalize_by_clusters(
    modes: list[Mode], x: np.ndarray, y: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """L2-orthonormalize raw modes inside each equal-eigenvalue cluster."""
    point_count = len(weights)
    mode_count = len(modes)
    values = np.empty((point_count, mode_count))
    gx = np.empty_like(values)
    gy = np.empty_like(values)
    for col, mode in enumerate(modes):
        values[:, col], gx[:, col], gy[:, col] = raw_lame_mode(mode, x, y)

    start = 0
    while start < mode_count:
        q = modes[start].q
        end = start + 1
        while end < mode_count and modes[end].q == q:
            end += 1
        block = slice(start, end)
        gram = values[:, block].T @ (weights[:, None] * values[:, block])
        transform = np.linalg.inv(np.linalg.cholesky(gram).T)
        values[:, block] = values[:, block] @ transform
        gx[:, block] = gx[:, block] @ transform
        gy[:, block] = gy[:, block] @ transform
        start = end
    return values, gx, gy


def complete_cluster_indices(modes: list[Mode], max_index: int) -> list[int]:
    """Return 1-based cutoffs that do not split an equal-eigenvalue cluster."""
    out: list[int] = []
    for n in range(1, max_index + 1):
        if n == len(modes) or modes[n - 1].q != modes[n].q:
            out.append(n)
    return out


def p_derivative_matrices() -> dict[str, np.ndarray]:
    """Derivatives of P(s,t)=S^{-1}S^{-T} at the equilateral point.

    Paper formula: P(s,t)=S^{-1}S^{-T} for the pullback to the equilateral triangle.
    Taylor derivatives of P are evaluated at p0.
    """
    y = Y0
    return {
        "x": np.array([[0.0, -1.0 / y], [-1.0 / y, 0.0]]),
        "y": np.array([[0.0, 0.0], [0.0, -2.0 / y]]),
        "xx": np.array([[2.0 / y**2, 0.0], [0.0, 0.0]]),
        "xy": np.array([[0.0, 2.0 / y**2], [2.0 / y**2, 0.0]]),
        "yy": np.array([[0.0, 0.0], [0.0, 6.0 / y**2]]),
        "xxx": np.zeros((2, 2)),
        "xxy": np.array([[-4.0 / y**3, 0.0], [0.0, 0.0]]),
    }


def mat_norm_2(mat: np.ndarray) -> float:
    return float(np.linalg.norm(mat, 2))


def apply_matrix_to_gradient(mat: np.ndarray, gx: np.ndarray, gy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return mat[0, 0] * gx + mat[0, 1] * gy, mat[1, 0] * gx + mat[1, 1] * gy


def gradient_pairing(
    fx: np.ndarray, fy: np.ndarray, gx: np.ndarray, gy: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    return (fx @ (weights[:, None] * gx)) + (fy @ (weights[:, None] * gy))


def build_spectral_probe(
    J: int,
    modes_d: list[Mode],
    gx_d: np.ndarray,
    gy_d: np.ndarray,
    weights: np.ndarray,
    tail_mode: str,
) -> SpectralProbe:
    """Compute lambda_xx(0,0), M_xxx, and M_xxy for one truncation level."""
    mats = p_derivative_matrices()
    lambdas = np.array([mode.eigenvalue for mode in modes_d])
    lambda1 = lambdas[0]
    grad1_norm = math.sqrt(lambda1)
    phi1_x = gx_d[:, 0]
    phi1_y = gy_d[:, 0]

    fx_x, fy_x = apply_matrix_to_gradient(mats["x"], phi1_x, phi1_y)
    fx_y, fy_y = apply_matrix_to_gradient(mats["y"], phi1_x, phi1_y)
    b_x = gradient_pairing(fx_x, fy_x, gx_d, gy_d, weights)
    b_y = gradient_pairing(fx_y, fy_y, gx_d, gy_d, weights)
    lambda_x = float(b_x[0])
    lambda_y = float(b_y[0])

    norm_fx_x = float(np.sum(weights * (fx_x**2 + fy_x**2)))
    norm_fx_y = float(np.sum(weights * (fx_y**2 + fy_y**2)))

    def w_bounds(b_vec: np.ndarray, forcing_norm_sq: float) -> tuple[float, float]:
        # Paper formula: w_a = -sum_{j>=2} B_j(a)/(lambda_j-lambda_1) * phi_j.
        # We bound ||w_a|| and ||grad w_a|| by finite spectral sums plus a tail bound.
        deltas = lambdas[1:J] - lambda1
        finite_b = b_vec[1:J]
        w0_sq = float(np.sum((finite_b**2) / (deltas**2)))
        w1_sq = float(np.sum(lambdas[1:J] * (finite_b**2) / (deltas**2)))
        if tail_mode == "dirichlet_parseval_tail":
            projected = float(np.sum((b_vec[:J] ** 2) / lambdas[:J]))
            tail_energy = max(forcing_norm_sq - projected, 0.0)
            lambda_next = lambdas[J]
            gap_next = lambda_next - lambda1
            w0_sq += lambda_next / (gap_next**2) * tail_energy
            w1_sq += (lambda_next**2) / (gap_next**2) * tail_energy
        elif tail_mode != "tail_ignored":
            raise ValueError(f"unknown tail mode: {tail_mode}")
        return math.sqrt(max(w0_sq, 0.0)), math.sqrt(max(w1_sq, 0.0))

    w0_x, w1_x = w_bounds(b_x, norm_fx_x)
    w0_y, w1_y = w_bounds(b_y, norm_fx_y)

    # Paper formula: spectral formula for lambda_xx at the equilateral triangle.
    axx_x, axx_y = apply_matrix_to_gradient(mats["xx"], phi1_x, phi1_y)
    a_xx = float(np.sum(weights * (axx_x * phi1_x + axx_y * phi1_y)))
    finite_second = float(np.sum((b_x[1:J] ** 2) / (lambdas[1:J] - lambda1)))
    lambda_xx0_lower = a_xx - 2.0 * finite_second
    if tail_mode == "dirichlet_parseval_tail":
        projected_x = float(np.sum((b_x[:J] ** 2) / lambdas[:J]))
        tail_x = max(norm_fx_x - projected_x, 0.0)
        lambda_next = lambdas[J]
        alpha = lambda_next / (lambda_next - lambda1)
        lambda_xx0_lower -= 2.0 * alpha * tail_x
    elif tail_mode != "tail_ignored":
        raise ValueError(f"unknown tail mode: {tail_mode}")

    def third_bound(d1: str, d2: str, d3: str) -> float:
        """Norm bound for the third derivative lambda_{d1,d2,d3}."""
        dirs = [d1, d2, d3]
        w0 = {"x": w0_x, "y": w0_y}
        w1 = {"x": w1_x, "y": w1_y}
        lam_bound = {
            "x": mat_norm_2(mats["x"]) * lambda1,
            "y": mat_norm_2(mats["y"]) * lambda1,
        }

        key3 = "".join(dirs)
        key3 = "".join(sorted(key3))
        if key3 == "xxx":
            p_abc = mats["xxx"]
        elif key3 == "xxy":
            p_abc = mats["xxy"]
        else:
            raise ValueError(f"third derivative key not implemented: {key3}")

        def second_key(a: str, b: str) -> str:
            return "".join(sorted(a + b))

        # Paper formula: third derivative without the second material derivative.
        # We use the norm upper bound to compute M_xxx and M_xxy.
        total = mat_norm_2(p_abc) * lambda1
        second_pairs = [
            (d1, d2, d3),
            (d1, d3, d2),
            (d2, d3, d1),
        ]
        for a, b, c in second_pairs:
            # The factor 4 follows the conservative first-experiment bound in the prompt.
            total += 4.0 * mat_norm_2(mats[second_key(a, b)]) * grad1_norm * w1[c]

        first_terms = [
            (d1, d2, d3),
            (d2, d1, d3),
            (d3, d1, d2),
        ]
        for a, b, c in first_terms:
            total += 2.0 * mat_norm_2(mats[a]) * w1[b] * w1[c]
            total += 2.0 * lam_bound[a] * w0[b] * w0[c]
        return float(total)

    return SpectralProbe(
        J=J,
        lambda1=lambda1,
        lambda_x=lambda_x,
        lambda_y=lambda_y,
        lambda_xx0_lower=float(lambda_xx0_lower),
        m_xxx=third_bound("x", "x", "x"),
        m_xxy=third_bound("x", "x", "y"),
        tail_mode=tail_mode,
        rigorous_flag="exploratory_double",
    )


def rxx_value(kind: str, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    r1 = np.sqrt(x**2 + y**2)
    r2 = np.sqrt((x - 1.0) ** 2 + y**2)
    perimeter_factor = 1.0 + r1 + r2
    base = (x / r1 + (x - 1.0) / r2) ** 2
    curvature = y**2 * (1.0 / r1**3 + 1.0 / r2**3)
    if kind == "J1":
        return -(math.pi**2 / (4.0 * y)) * (base + perimeter_factor * curvature)
    if kind == "J2":
        c_star = 4.0 * math.pi**2 / (3.0 + math.sqrt(math.pi * math.sqrt(3.0))) ** 2
        q = perimeter_factor + np.sqrt(2.0 * math.pi * y)
        return -(c_star / y) * (base + q * curvature)
    raise ValueError(f"unknown functional: {kind}")


def rxx_lower_sampled(kind: str, rx: float, ry: float, samples: int) -> float:
    """Exploratory sampled lower bound for Rxx over the cell."""
    ss = np.linspace(0.0, rx, samples)
    tt = np.linspace(-ry, 0.0, samples)
    s_grid, t_grid = np.meshgrid(ss, tt, indexing="ij")
    x = 0.5 + s_grid
    y = Y0 + t_grid
    return float(np.min(rxx_value(kind, x, y)))


def jxx_lower(probe: SpectralProbe, rx: float, ry: float, samples: int) -> tuple[float, float, float, float]:
    """Evaluate the Taylor lower bound for both J1_xx and J2_xx.

    Paper formula: Taylor lower bound for J_xx near the equilateral triangle.
    The cell is accepted if both L_1 and L_2 are positive.
    """
    lambda_xx_cell_lower = probe.lambda_xx0_lower - rx * probe.m_xxx - ry * probe.m_xxy
    spectral_part = 0.5 * (Y0 - ry) * lambda_xx_cell_lower
    rxx1 = rxx_lower_sampled("J1", rx, ry, samples)
    rxx2 = rxx_lower_sampled("J2", rx, ry, samples)
    return spectral_part + rxx1, spectral_part + rxx2, rxx1, rxx2


def shape_radii(cell_shape: str, radius: float) -> tuple[float, float]:
    if cell_shape == "x_only":
        return radius, 0.0
    if cell_shape == "y_only":
        return 0.0, radius
    if cell_shape == "omega_up_type":
        return radius, 0.1 * radius
    if cell_shape == "square":
        return radius, radius
    raise ValueError(f"unknown cell shape: {cell_shape}")


def accepted(probe: SpectralProbe, cell_shape: str, radius: float, samples: int) -> tuple[bool, tuple[float, ...]]:
    rx, ry = shape_radii(cell_shape, radius)
    if ry >= Y0:
        return False, (float("-inf"), float("-inf"), float("nan"), float("nan"))
    j1, j2, r1, r2 = jxx_lower(probe, rx, ry, samples)
    return j1 > 0.0 and j2 > 0.0, (j1, j2, r1, r2)


def bisect_cell_size(probe: SpectralProbe, cell_shape: str, max_radius: float, samples: int) -> tuple[float, tuple[float, ...]]:
    ok0, _ = accepted(probe, cell_shape, 0.0, samples)
    if not ok0:
        return 0.0, jxx_lower(probe, *shape_radii(cell_shape, 0.0), samples)

    lo = 0.0
    hi = max_radius
    ok_hi, _ = accepted(probe, cell_shape, hi, samples)
    if ok_hi:
        lo = hi
    else:
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            ok_mid, _ = accepted(probe, cell_shape, mid, samples)
            if ok_mid:
                lo = mid
            else:
                hi = mid
    rx, ry = shape_radii(cell_shape, lo)
    return lo, jxx_lower(probe, rx, ry, samples)


def write_outputs(
    outdir: Path,
    rows: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    max_csv = outdir / "max_cell_sizes.csv"
    diag_csv = outdir / "diagnostics.csv"
    fields = [
        "model",
        "J",
        "tail_mode",
        "cell_shape",
        "rx",
        "ry",
        "distance",
        "lambda_xx0",
        "M_xxx",
        "M_xxy",
        "Rxx1_lower",
        "Rxx2_lower",
        "Jxx1_lower",
        "Jxx2_lower",
        "success",
        "rigorous_flag",
    ]
    with max_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with diag_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(diagnostics)

    summary = outdir / "summary.md"
    with summary.open("w") as f:
        f.write("# Omega Up Taylor Probe\n\n")
        f.write("This is an exploratory double-precision Taylor spectral probe. It does not modify the production Omega_up or Omega_mid verification code.\n\n")
        f.write("The spectral data use analytic equilateral triangle eigenfunctions, but spectral inner products and Rxx cell lower bounds are evaluated numerically in double precision. Therefore the rigorous flag is `exploratory_double`.\n\n")
        f.write("Known comparison points:\n\n")
        f.write("- Current cellwise eigenspace enclosure succeeds only around distance `1e-5`.\n")
        f.write("- Idealized fixed-equilateral-coefficient model suggests `1e-2` to `1e-1` scale.\n")
        f.write("- This Taylor probe tests whether the third-derivative Taylor certificate can bridge that gap.\n\n")
        f.write("| model | J | tail mode | cell shape | max rx | max ry | distance | Jxx1 lower | Jxx2 lower | rigorous flag |\n")
        f.write("|---|---:|---|---|---:|---:|---:|---:|---:|---|\n")
        for row in rows:
            f.write(
                f"| {row['model']} | {row['J']} | {row['tail_mode']} | {row['cell_shape']} | "
                f"{float(row['rx']):.8e} | {float(row['ry']):.8e} | {float(row['distance']):.8e} | "
                f"{float(row['Jxx1_lower']):.8e} | {float(row['Jxx2_lower']):.8e} | {row['rigorous_flag']} |\n"
            )

        omega_rows = [row for row in rows if row["cell_shape"] == "omega_up_type"]
        if omega_rows:
            best = max(omega_rows, key=lambda r: float(r["rx"]))
            f.write("\nMost important Omega_up-type result:\n\n")
            f.write(
                f"- best `rx = {float(best['rx']):.8e}`, `ry = {float(best['ry']):.8e}` "
                f"at `J = {best['J']}`, tail mode `{best['tail_mode']}`.\n"
            )
            f.write(
                f"- target `rx >= 1e-2`, `ry = 1e-3`: "
                f"{'reached' if float(best['rx']) >= 1e-2 and float(best['ry']) >= 1e-3 else 'not reached'}.\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--J-list", default="30,50,100,200", help="comma-separated truncation levels")
    parser.add_argument(
        "--tail-modes",
        default="tail_ignored,dirichlet_parseval_tail",
        help="comma-separated tail modes: tail_ignored, dirichlet_parseval_tail",
    )
    parser.add_argument("--nq", type=int, default=260, help="1D Gauss quadrature order")
    parser.add_argument("--rxx-samples", type=int, default=81, help="sample count per direction for Rxx lower diagnostic")
    parser.add_argument("--max-radius", type=float, default=0.2, help="largest radius searched in bisection")
    parser.add_argument("--outdir", default="results/omega_up_taylor_probe", help="output directory")
    args = parser.parse_args()

    j_list = parse_j_list(args.J_list)
    tail_modes = [x.strip() for x in args.tail_modes.split(",") if x.strip()]
    max_j = max(j_list)
    needed_modes = max_j + 2
    x, y, weights = triangle_quadrature(args.nq)
    modes_d = enumerate_equilateral_modes("D", needed_modes)
    _, gx_d, gy_d = orthonormalize_by_clusters(modes_d, x, y, weights)

    rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    cell_shapes = ["x_only", "y_only", "omega_up_type", "square"]
    target_checks = [
        ("omega_up_target_1e-2", 1e-2, 1e-3),
    ]

    for j in j_list:
        for tail_mode in tail_modes:
            probe = build_spectral_probe(j, modes_d, gx_d, gy_d, weights, tail_mode)
            print(
                f"J={j}, tail_mode={tail_mode}: lambda_xx0_lower={probe.lambda_xx0_lower:.8e}, "
                f"M_xxx={probe.m_xxx:.8e}, M_xxy={probe.m_xxy:.8e}",
                flush=True,
            )
            for shape in cell_shapes:
                radius, values = bisect_cell_size(probe, shape, args.max_radius, args.rxx_samples)
                rx, ry = shape_radii(shape, radius)
                j1, j2, r1, r2 = values
                rows.append(
                    {
                        "model": "taylor_spectral_probe",
                        "J": j,
                        "tail_mode": probe.tail_mode,
                        "cell_shape": shape,
                        "rx": rx,
                        "ry": ry,
                        "distance": math.hypot(rx, ry),
                        "lambda_xx0": probe.lambda_xx0_lower,
                        "M_xxx": probe.m_xxx,
                        "M_xxy": probe.m_xxy,
                        "Rxx1_lower": r1,
                        "Rxx2_lower": r2,
                        "Jxx1_lower": j1,
                        "Jxx2_lower": j2,
                        "success": j1 > 0.0 and j2 > 0.0,
                        "rigorous_flag": probe.rigorous_flag,
                    }
                )

            for label, rx, ry in target_checks:
                j1, j2, r1, r2 = jxx_lower(probe, rx, ry, args.rxx_samples)
                diagnostics.append(
                {
                    "model": "taylor_spectral_probe",
                    "J": j,
                    "tail_mode": probe.tail_mode,
                        "cell_shape": label,
                        "rx": rx,
                        "ry": ry,
                        "distance": math.hypot(rx, ry),
                        "lambda_xx0": probe.lambda_xx0_lower,
                        "M_xxx": probe.m_xxx,
                        "M_xxy": probe.m_xxy,
                        "Rxx1_lower": r1,
                        "Rxx2_lower": r2,
                        "Jxx1_lower": j1,
                        "Jxx2_lower": j2,
                        "success": j1 > 0.0 and j2 > 0.0,
                        "rigorous_flag": probe.rigorous_flag,
                    }
                )

    write_outputs(Path(args.outdir), rows, diagnostics)
    print(f"Wrote {args.outdir}/summary.md")
    print(f"Wrote {args.outdir}/max_cell_sizes.csv")
    print(f"Wrote {args.outdir}/diagnostics.csv")


if __name__ == "__main__":
    main()
