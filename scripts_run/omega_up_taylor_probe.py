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
    lambda_xx_finite: float
    lambda_xxx_finite: float
    lambda_xxy_finite: float
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

    # Paper formula: finite spectral formula for lambda_xx at the equilateral triangle.
    axx_x, axx_y = apply_matrix_to_gradient(mats["xx"], phi1_x, phi1_y)
    a_xx = float(np.sum(weights * (axx_x * phi1_x + axx_y * phi1_y)))
    finite_second = float(np.sum((b_x[1:J] ** 2) / (lambdas[1:J] - lambda1)))
    lambda_xx_finite = a_xx - 2.0 * finite_second
    lambda_xx0_lower = lambda_xx_finite
    if tail_mode == "dirichlet_parseval_tail":
        projected_x = float(np.sum((b_x[:J] ** 2) / lambdas[:J]))
        tail_x = max(norm_fx_x - projected_x, 0.0)
        lambda_next = lambdas[J]
        alpha = lambda_next / (lambda_next - lambda1)
        lambda_xx0_lower -= 2.0 * alpha * tail_x
    elif tail_mode != "tail_ignored":
        raise ValueError(f"unknown tail mode: {tail_mode}")

    def material_coefficients(b_vec: np.ndarray) -> np.ndarray:
        coeff = np.zeros(J)
        coeff[1:J] = -b_vec[1:J] / (lambdas[1:J] - lambda1)
        return coeff

    coeff_x = material_coefficients(b_x)
    coeff_y = material_coefficients(b_y)
    wx_gx = gx_d[:, :J] @ coeff_x
    wx_gy = gy_d[:, :J] @ coeff_x
    wy_gx = gx_d[:, :J] @ coeff_y
    wy_gy = gy_d[:, :J] @ coeff_y

    def second_key(a: str, b: str) -> str:
        return "".join(sorted(a + b))

    def third_key(a: str, b: str, c: str) -> str:
        return "".join(sorted(a + b + c))

    def third_matrix(key: str) -> np.ndarray:
        if key not in {"xxx", "xxy"}:
            raise ValueError(f"third derivative key not implemented: {key}")
        return mats[key]

    def grad_inner_matrix(
        mat: np.ndarray,
        ax: np.ndarray,
        ay: np.ndarray,
        bx: np.ndarray,
        by: np.ndarray,
    ) -> float:
        mx, my = apply_matrix_to_gradient(mat, ax, ay)
        return float(np.sum(weights * (mx * bx + my * by)))

    def l2_inner_coeffs(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.dot(a, b))

    def third_finite(d1: str, d2: str, d3: str) -> float:
        """Finite signed third derivative using the paper formula without u_ab."""
        grad_w = {
            "x": (wx_gx, wx_gy),
            "y": (wy_gx, wy_gy),
        }
        coeff_w = {
            "x": coeff_x,
            "y": coeff_y,
        }
        lambda_dir = {
            "x": lambda_x,
            "y": lambda_y,
        }

        total = grad_inner_matrix(third_matrix(third_key(d1, d2, d3)), phi1_x, phi1_y, phi1_x, phi1_y)
        second_pairs = [
            (d1, d2, d3),
            (d1, d3, d2),
            (d2, d3, d1),
        ]
        for a, b, c in second_pairs:
            c_gx, c_gy = grad_w[c]
            total += 2.0 * grad_inner_matrix(mats[second_key(a, b)], c_gx, c_gy, phi1_x, phi1_y)

        first_terms = [
            (d1, d2, d3),
            (d2, d1, d3),
            (d3, d1, d2),
        ]
        for a, b, c in first_terms:
            b_gx, b_gy = grad_w[b]
            c_gx, c_gy = grad_w[c]
            total += 2.0 * grad_inner_matrix(mats[a], b_gx, b_gy, c_gx, c_gy)
            total -= 2.0 * lambda_dir[a] * l2_inner_coeffs(coeff_w[b], coeff_w[c])
        return float(total)

    def third_bound(d1: str, d2: str, d3: str) -> float:
        """Norm bound for the third derivative lambda_{d1,d2,d3}."""
        dirs = [d1, d2, d3]
        w0 = {"x": w0_x, "y": w0_y}
        w1 = {"x": w1_x, "y": w1_y}
        lam_bound = {
            "x": mat_norm_2(mats["x"]) * lambda1,
            "y": mat_norm_2(mats["y"]) * lambda1,
        }

        p_abc = third_matrix(third_key(*dirs))

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
        lambda_xx_finite=float(lambda_xx_finite),
        lambda_xxx_finite=third_finite("x", "x", "x"),
        lambda_xxy_finite=third_finite("x", "x", "y"),
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


def lambda_xx_cell_lower(probe: SpectralProbe, rx: float, ry: float, model: str) -> float:
    """Lower model for lambda_xx over the near-equilateral cell."""
    if model == "norm_bound":
        return probe.lambda_xx0_lower - rx * probe.m_xxx - ry * probe.m_xxy
    if model == "symmetry_first_order":
        return probe.lambda_xx_finite - ry * probe.m_xxy
    if model == "signed_xxy_first_order":
        return probe.lambda_xx_finite + min(-probe.lambda_xxy_finite * ry, 0.0)
    raise ValueError(f"unknown Taylor model: {model}")


def jxx_lower(
    probe: SpectralProbe,
    rx: float,
    ry: float,
    samples: int,
    model: str,
) -> tuple[float, float, float, float, float]:
    """Evaluate the Taylor lower bound for both J1_xx and J2_xx.

    Paper formula: Taylor lower bound for J_xx near the equilateral triangle.
    The cell is accepted if both L_1 and L_2 are positive.
    """
    cell_lambda_xx = lambda_xx_cell_lower(probe, rx, ry, model)
    spectral_part = 0.5 * (Y0 - ry) * cell_lambda_xx
    rxx1 = rxx_lower_sampled("J1", rx, ry, samples)
    rxx2 = rxx_lower_sampled("J2", rx, ry, samples)
    return spectral_part + rxx1, spectral_part + rxx2, rxx1, rxx2, cell_lambda_xx


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


def accepted(
    probe: SpectralProbe,
    cell_shape: str,
    radius: float,
    samples: int,
    model: str,
) -> tuple[bool, tuple[float, ...]]:
    rx, ry = shape_radii(cell_shape, radius)
    if ry >= Y0:
        return False, (float("-inf"), float("-inf"), float("nan"), float("nan"), float("nan"))
    j1, j2, r1, r2, cell_lambda_xx = jxx_lower(probe, rx, ry, samples, model)
    return j1 > 0.0 and j2 > 0.0, (j1, j2, r1, r2, cell_lambda_xx)


def bisect_cell_size(
    probe: SpectralProbe,
    cell_shape: str,
    max_radius: float,
    samples: int,
    model: str,
) -> tuple[float, tuple[float, ...]]:
    ok0, _ = accepted(probe, cell_shape, 0.0, samples, model)
    if not ok0:
        return 0.0, jxx_lower(probe, *shape_radii(cell_shape, 0.0), samples, model)

    lo = 0.0
    hi = max_radius
    ok_hi, _ = accepted(probe, cell_shape, hi, samples, model)
    if ok_hi:
        lo = hi
    else:
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            ok_mid, _ = accepted(probe, cell_shape, mid, samples, model)
            if ok_mid:
                lo = mid
            else:
                hi = mid
    rx, ry = shape_radii(cell_shape, lo)
    return lo, jxx_lower(probe, rx, ry, samples, model)


def cluster_complete(modes: list[Mode], J: int) -> bool:
    return J == len(modes) or modes[J - 1].q != modes[J].q


def neighboring_complete_cutoffs(modes: list[Mode], J: int) -> tuple[int, int]:
    """Nearest complete 1-based cluster cutoffs around the requested cutoff J."""
    cutoffs = complete_cluster_indices(modes, len(modes))
    prev_complete = max(c for c in cutoffs if c <= J)
    next_complete = min(c for c in cutoffs if c >= J)
    return prev_complete, next_complete


def lambda_xx_shell_contribution_rows(
    modes_d: list[Mode],
    gx_d: np.ndarray,
    gy_d: np.ndarray,
    weights: np.ndarray,
    max_j: int,
) -> list[dict[str, object]]:
    """Shell-wise finite sum for lambda_xx at the equilateral triangle.

    Paper formula:
        lambda_xx = (P_xx grad phi_1, grad phi_1)
                    - 2 sum_{j>=2} B_j^2/(lambda_j-lambda_1).

    Equal-eigenvalue shells are summed together because individual basis
    coefficients inside a repeated eigenspace are basis-dependent.
    """
    mats = p_derivative_matrices()
    lambdas = np.array([mode.eigenvalue for mode in modes_d])
    lambda1 = lambdas[0]
    phi1_x = gx_d[:, 0]
    phi1_y = gy_d[:, 0]
    g_x, g_y = apply_matrix_to_gradient(mats["x"], phi1_x, phi1_y)
    b_x = gradient_pairing(g_x, g_y, gx_d, gy_d, weights)
    axx_x, axx_y = apply_matrix_to_gradient(mats["xx"], phi1_x, phi1_y)
    a_xx = float(np.sum(weights * (axx_x * phi1_x + axx_y * phi1_y)))

    rows: list[dict[str, object]] = []
    cumulative = 0.0
    cumulative_parseval = 0.0
    shell_index = 0
    start = 1
    while start < len(modes_d):
        q = modes_d[start].q
        end = start + 1
        while end < len(modes_d) and modes_d[end].q == q:
            end += 1
        if start >= max_j:
            break
        shell_index += 1
        shell_lambda = lambdas[start]
        shell_b_sq = float(np.sum(b_x[start:end] ** 2))
        contribution = shell_b_sq / (shell_lambda - lambda1)
        parseval_contribution = shell_b_sq / shell_lambda
        cumulative += contribution
        cumulative_parseval += parseval_contribution
        rows.append(
            {
                "shell_index": shell_index,
                "q": q,
                "eigenvalue": shell_lambda,
                "multiplicity": end - start,
                "mode_start": start + 1,
                "mode_end": end,
                "complete_J": end,
                "requested_max_J_cuts_shell": start < max_j < end,
                "B_sq_shell": shell_b_sq,
                "contribution": contribution,
                "cumulative": cumulative,
                "lambda_xx_finite_shell": a_xx - 2.0 * cumulative,
                "parseval_contribution": parseval_contribution,
                "parseval_cumulative": cumulative_parseval,
                "rigorous_flag": "exploratory_double",
            }
        )
        start = end

    final_cumulative = rows[-1]["cumulative"] if rows else 0.0
    for row in rows:
        row["fraction_of_last_cumulative"] = float(row["cumulative"]) / final_cumulative if final_cumulative else float("nan")
    return rows


def basis_orthogonality_rows(
    j_list: list[int],
    bc: str,
    modes: list[Mode],
    values: np.ndarray,
    gx: np.ndarray,
    gy: np.ndarray,
    weights: np.ndarray,
) -> list[dict[str, object]]:
    """Check L2 orthonormality and stiffness diagonality for analytic modes."""
    max_j = max(j_list)
    lambdas = np.array([mode.eigenvalue for mode in modes[:max_j]])
    values_j = values[:, :max_j]
    gx_j = gx[:, :max_j]
    gy_j = gy[:, :max_j]
    mass = values_j.T @ (weights[:, None] * values_j)
    stiffness = gx_j.T @ (weights[:, None] * gx_j) + gy_j.T @ (weights[:, None] * gy_j)
    eye = np.eye(max_j)
    lambda_diag = np.diag(lambdas)
    offdiag_mask = ~np.eye(max_j, dtype=bool)

    rows: list[dict[str, object]] = []
    for J in j_list:
        mass_j = mass[:J, :J]
        stiffness_j = stiffness[:J, :J]
        mass_error = mass_j - eye[:J, :J]
        stiffness_error = stiffness_j - lambda_diag[:J, :J]
        offdiag_j = offdiag_mask[:J, :J]
        rows.append(
            {
                "bc": bc,
                "J": J,
                "cluster_complete": cluster_complete(modes, J),
                "mass_max_abs_error": float(np.max(np.abs(mass_error))),
                "mass_diag_max_abs_error": float(np.max(np.abs(np.diag(mass_j) - 1.0))),
                "mass_offdiag_max_abs": float(np.max(np.abs(mass_j[offdiag_j]))) if J > 1 else 0.0,
                "stiffness_max_abs_error": float(np.max(np.abs(stiffness_error))),
                "stiffness_relative_max_error": float(np.max(np.abs(stiffness_error)) / np.max(lambdas[:J])),
                "stiffness_diag_relative_max_error": float(
                    np.max(np.abs((np.diag(stiffness_j) - lambdas[:J]) / lambdas[:J]))
                ),
                "stiffness_offdiag_max_abs": float(np.max(np.abs(stiffness_j[offdiag_j]))) if J > 1 else 0.0,
                "rigorous_flag": "exploratory_double",
            }
        )
    return rows


def parseval_decomposition_rows(
    j_list: list[int],
    modes_d: list[Mode],
    gx_d: np.ndarray,
    gy_d: np.ndarray,
    modes_n: list[Mode],
    gx_n: np.ndarray,
    gy_n: np.ndarray,
    weights: np.ndarray,
) -> list[dict[str, object]]:
    """Diagnostic Helmholtz/Parseval split of G=P_x grad(phi_1)."""
    mats = p_derivative_matrices()
    lambda_d = np.array([mode.eigenvalue for mode in modes_d])
    mu_n = np.array([mode.eigenvalue for mode in modes_n])
    g_x, g_y = apply_matrix_to_gradient(mats["x"], gx_d[:, 0], gy_d[:, 0])
    g_norm_sq = float(np.sum(weights * (g_x**2 + g_y**2)))

    dirichlet_coupling = gradient_pairing(g_x, g_y, gx_d, gy_d, weights)
    # Paper diagnostic: Neumann rotated-gradient part, (G, J grad psi_m).
    neumann_rotated_coupling = (g_x @ (weights[:, None] * (-gy_n))) + (g_y @ (weights[:, None] * gx_n))

    d_sums = {
        J: float(np.sum((dirichlet_coupling[:J] ** 2) / lambda_d[:J]))
        for J in j_list
    }
    n_sums = {
        M: float(np.sum((neumann_rotated_coupling[:M] ** 2) / mu_n[:M]))
        for M in j_list
    }

    rows: list[dict[str, object]] = []
    for J in j_list:
        for M in j_list:
            d_j = d_sums[J]
            n_m = n_sums[M]
            combined = d_j + n_m
            rows.append(
                {
                    "J": J,
                    "M": M,
                    "D_J": d_j,
                    "N_M": n_m,
                    "G_norm_sq": g_norm_sq,
                    "dirichlet_residual": g_norm_sq - d_j,
                    "combined_residual": g_norm_sq - combined,
                    "D_fraction": d_j / g_norm_sq,
                    "N_fraction": n_m / g_norm_sq,
                    "combined_fraction": combined / g_norm_sq,
                    "rigorous_flag": "exploratory_double",
                }
            )
    return rows


def write_outputs(
    outdir: Path,
    rows: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
    third_rows: list[dict[str, object]],
    parseval_rows: list[dict[str, object]],
    shell_rows: list[dict[str, object]],
    orthogonality_rows: list[dict[str, object]],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    max_csv = outdir / "max_cell_sizes.csv"
    diag_csv = outdir / "diagnostics.csv"
    third_csv = outdir / "third_derivative_diagnostics.csv"
    parseval_csv = outdir / "parseval_decomposition.csv"
    shell_csv = outdir / "lambda_xx_shell_contributions.csv"
    orthogonality_csv = outdir / "basis_orthogonality_diagnostics.csv"
    fields = [
        "model",
        "J",
        "tail_mode",
        "cell_shape",
        "rx",
        "ry",
        "distance",
        "lambda_xx_cell_lower",
        "lambda_xx0",
        "lambda_xx_finite",
        "lambda_xxx_finite",
        "lambda_xxy_finite",
        "M_xxx",
        "M_xxy",
        "Rxx1_lower",
        "Rxx2_lower",
        "Jxx1_lower",
        "Jxx2_lower",
        "success",
        "cluster_complete",
        "rigorous_flag",
    ]
    third_fields = [
        "J",
        "cluster_complete",
        "prev_complete_J",
        "next_complete_J",
        "lambda_x_finite",
        "lambda_y_finite",
        "lambda_xx_finite",
        "lambda_xxx_finite",
        "lambda_xxy_finite",
        "abs_lambda_xxx_finite",
        "M_xxx_norm",
        "ratio_xxx",
        "abs_lambda_xxy_finite",
        "M_xxy_norm",
        "ratio_xxy",
        "rigorous_flag",
    ]
    parseval_fields = [
        "J",
        "M",
        "D_J",
        "N_M",
        "G_norm_sq",
        "dirichlet_residual",
        "combined_residual",
        "D_fraction",
        "N_fraction",
        "combined_fraction",
        "rigorous_flag",
    ]
    shell_fields = [
        "shell_index",
        "q",
        "eigenvalue",
        "multiplicity",
        "mode_start",
        "mode_end",
        "complete_J",
        "requested_max_J_cuts_shell",
        "B_sq_shell",
        "contribution",
        "cumulative",
        "fraction_of_last_cumulative",
        "lambda_xx_finite_shell",
        "parseval_contribution",
        "parseval_cumulative",
        "rigorous_flag",
    ]
    orthogonality_fields = [
        "bc",
        "J",
        "cluster_complete",
        "mass_max_abs_error",
        "mass_diag_max_abs_error",
        "mass_offdiag_max_abs",
        "stiffness_max_abs_error",
        "stiffness_relative_max_error",
        "stiffness_diag_relative_max_error",
        "stiffness_offdiag_max_abs",
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
    with third_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=third_fields)
        writer.writeheader()
        writer.writerows(third_rows)
    with parseval_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parseval_fields)
        writer.writeheader()
        writer.writerows(parseval_rows)
    with shell_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=shell_fields)
        writer.writeheader()
        writer.writerows(shell_rows)
    with orthogonality_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=orthogonality_fields)
        writer.writeheader()
        writer.writerows(orthogonality_rows)

    summary = outdir / "summary.md"
    with summary.open("w") as f:
        f.write("# Omega Up Taylor Probe\n\n")
        f.write("This is an exploratory double-precision Taylor spectral probe. It does not modify the production Omega_up or Omega_mid verification code.\n\n")
        f.write("The spectral data use analytic equilateral triangle eigenfunctions, but spectral inner products and Rxx cell lower bounds are evaluated numerically in double precision. Therefore the rigorous flag is `exploratory_double`.\n\n")
        f.write("Known comparison points:\n\n")
        f.write("- Current cellwise eigenspace enclosure succeeds only around distance `1e-5`.\n")
        f.write("- Idealized fixed-equilateral-coefficient model suggests `1e-2` to `1e-1` scale.\n")
        f.write("- This Taylor probe tests whether the third-derivative Taylor certificate can bridge that gap.\n\n")
        f.write("A row that reaches the configured `--max-radius` search cap should be read as a lower diagnostic, not as the true maximal accepted radius.\n\n")
        f.write("Taylor models:\n\n")
        f.write("- `norm_bound`: existing first-order Taylor lower bound using `-rx M_xxx - ry M_xxy`.\n")
        f.write("- `symmetry_first_order`: finite `lambda_xx(0,0)` with the `rx M_xxx` penalty removed by equilateral x-symmetry, keeping `-ry M_xxy`.\n")
        f.write("- `signed_xxy_first_order`: finite `lambda_xx(0,0)` plus the signed affine `lambda_xxy t` minimum on `t in [-ry,0]`.\n")
        f.write("- `dirichlet_parseval_tail`: the intentionally crude Dirichlet-only Parseval tail diagnostic; it is not used by the two symmetry-improved models.\n\n")

        f.write("Finite signed derivative diagnostics:\n\n")
        f.write("| J | cluster complete | complete cutoffs | lambda_x | lambda_xx finite | lambda_xxx finite | lambda_xxy finite | M_xxx/abs | M_xxy/abs |\n")
        f.write("|---:|---|---|---:|---:|---:|---:|---:|---:|\n")
        for row in third_rows:
            f.write(
                f"| {row['J']} | {row['cluster_complete']} | "
                f"{row['prev_complete_J']}/{row['next_complete_J']} | "
                f"{float(row['lambda_x_finite']):.8e} | {float(row['lambda_xx_finite']):.8e} | "
                f"{float(row['lambda_xxx_finite']):.8e} | "
                f"{float(row['lambda_xxy_finite']):.8e} | {float(row['ratio_xxx']):.8e} | "
                f"{float(row['ratio_xxy']):.8e} |\n"
            )
        bad_xxx = [row for row in third_rows if abs(float(row["lambda_xxx_finite"])) > 1.0e-6]
        if bad_xxx:
            f.write("\nWarning: `lambda_xxx_finite` is not close to zero for at least one requested J. This can happen if the requested truncation cuts an equal-eigenvalue cluster, but it should be treated as an implementation/integration diagnostic.\n")
        else:
            f.write("\nThe finite `lambda_xxx` values are numerically close to zero at the `1e-6` level.\n")

        if orthogonality_rows:
            f.write("\nBasis orthogonality diagnostics:\n\n")
            f.write("| bc | J | mass max error | stiffness max error | stiffness relative error |\n")
            f.write("|---|---:|---:|---:|---:|\n")
            for row in orthogonality_rows:
                f.write(
                    f"| {row['bc']} | {row['J']} | {float(row['mass_max_abs_error']):.8e} | "
                    f"{float(row['stiffness_max_abs_error']):.8e} | "
                    f"{float(row['stiffness_relative_max_error']):.8e} |\n"
                )

        if shell_rows:
            f.write("\nShell-wise `lambda_xx` finite sum diagnostic:\n\n")
            f.write("| shell | q | multiplicity | complete J | contribution | cumulative | lambda_xx finite | fraction of final cumulative |\n")
            f.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in shell_rows[:12]:
                f.write(
                    f"| {row['shell_index']} | {row['q']} | {row['multiplicity']} | {row['complete_J']} | "
                    f"{float(row['contribution']):.8e} | {float(row['cumulative']):.8e} | "
                    f"{float(row['lambda_xx_finite_shell']):.8e} | "
                    f"{float(row['fraction_of_last_cumulative']):.8e} |\n"
                )
            f.write("- Full shell table is written to `lambda_xx_shell_contributions.csv`.\n")
        f.write("\nCell-size bisection results:\n\n")
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
                f"with model `{best['model']}`, `J = {best['J']}`, tail mode `{best['tail_mode']}`.\n"
            )
            f.write(
                f"- target `rx >= 1e-2`, `ry = 1e-3`: "
                f"{'reached' if float(best['rx']) >= 1e-2 and float(best['ry']) >= 1e-3 else 'not reached'}.\n"
            )
            f.write("\nOmega_up-type best row by model:\n\n")
            f.write("| model | best rx | best ry | J | tail mode | Jxx2 lower |\n")
            f.write("|---|---:|---:|---:|---|---:|\n")
            for model in sorted({row["model"] for row in omega_rows}):
                model_rows = [row for row in omega_rows if row["model"] == model]
                model_best = max(model_rows, key=lambda r: float(r["rx"]))
                f.write(
                    f"| {model} | {float(model_best['rx']):.8e} | {float(model_best['ry']):.8e} | "
                    f"{model_best['J']} | {model_best['tail_mode']} | {float(model_best['Jxx2_lower']):.8e} |\n"
                )

        target_rows = [row for row in diagnostics if row["cell_shape"] == "omega_up_target_1e-2"]
        if target_rows:
            f.write("\nTarget diagnostic at `rx=1e-2`, `ry=1e-3`:\n\n")
            f.write("| model | J | tail mode | Jxx1 lower | Jxx2 lower | success |\n")
            f.write("|---|---:|---|---:|---:|---|\n")
            for row in target_rows:
                f.write(
                    f"| {row['model']} | {row['J']} | {row['tail_mode']} | "
                    f"{float(row['Jxx1_lower']):.8e} | {float(row['Jxx2_lower']):.8e} | {row['success']} |\n"
                )

        if parseval_rows:
            best_parseval = min(parseval_rows, key=lambda r: abs(float(r["combined_residual"])))
            f.write("\nParseval decomposition diagnostic:\n\n")
            f.write(
                f"- smallest absolute combined residual in the grid is `{float(best_parseval['combined_residual']):.8e}` "
                f"at `J={best_parseval['J']}`, `M={best_parseval['M']}`.\n"
            )
            f.write("- Full values are written to `parseval_decomposition.csv`.\n")


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
    needed_modes = max_j + 128
    x, y, weights = triangle_quadrature(args.nq)
    modes_d = enumerate_equilateral_modes("D", needed_modes)
    values_d, gx_d, gy_d = orthonormalize_by_clusters(modes_d, x, y, weights)
    modes_n = enumerate_equilateral_modes("N", max_j + 128)
    values_n, gx_n, gy_n = orthonormalize_by_clusters(modes_n, x, y, weights)

    rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    third_rows: list[dict[str, object]] = []
    cell_shapes = ["x_only", "y_only", "omega_up_type", "square"]
    target_checks = [
        ("omega_up_target_1e-2", 1e-2, 1e-3),
    ]

    def output_row(
        model: str,
        probe: SpectralProbe,
        cell_shape: str,
        rx: float,
        ry: float,
        values: tuple[float, ...],
    ) -> dict[str, object]:
        j1, j2, r1, r2, cell_lambda_xx = values
        return {
            "model": model,
            "J": probe.J,
            "tail_mode": probe.tail_mode,
            "cell_shape": cell_shape,
            "rx": rx,
            "ry": ry,
            "distance": math.hypot(rx, ry),
            "lambda_xx_cell_lower": cell_lambda_xx,
            "lambda_xx0": probe.lambda_xx0_lower,
            "lambda_xx_finite": probe.lambda_xx_finite,
            "lambda_xxx_finite": probe.lambda_xxx_finite,
            "lambda_xxy_finite": probe.lambda_xxy_finite,
            "M_xxx": probe.m_xxx,
            "M_xxy": probe.m_xxy,
            "Rxx1_lower": r1,
            "Rxx2_lower": r2,
            "Jxx1_lower": j1,
            "Jxx2_lower": j2,
            "success": j1 > 0.0 and j2 > 0.0,
            "cluster_complete": cluster_complete(modes_d, probe.J),
            "rigorous_flag": probe.rigorous_flag,
        }

    for j in j_list:
        probe_tail_modes = list(dict.fromkeys(tail_modes + ["tail_ignored"]))
        probes: dict[str, SpectralProbe] = {}
        for tail_mode in probe_tail_modes:
            probe = build_spectral_probe(j, modes_d, gx_d, gy_d, weights, tail_mode)
            probes[tail_mode] = probe
            print(
                f"J={j}, tail_mode={tail_mode}: lambda_xx0_lower={probe.lambda_xx0_lower:.8e}, "
                f"lambda_xx_finite={probe.lambda_xx_finite:.8e}, "
                f"lambda_xxx_finite={probe.lambda_xxx_finite:.8e}, "
                f"lambda_xxy_finite={probe.lambda_xxy_finite:.8e}, "
                f"M_xxx={probe.m_xxx:.8e}, M_xxy={probe.m_xxy:.8e}",
                flush=True,
            )

        finite_probe = probes["tail_ignored"]
        prev_complete, next_complete = neighboring_complete_cutoffs(modes_d, j)
        third_rows.append(
            {
                "J": j,
                "cluster_complete": cluster_complete(modes_d, j),
                "prev_complete_J": prev_complete,
                "next_complete_J": next_complete,
                "lambda_x_finite": finite_probe.lambda_x,
                "lambda_y_finite": finite_probe.lambda_y,
                "lambda_xx_finite": finite_probe.lambda_xx_finite,
                "lambda_xxx_finite": finite_probe.lambda_xxx_finite,
                "lambda_xxy_finite": finite_probe.lambda_xxy_finite,
                "abs_lambda_xxx_finite": abs(finite_probe.lambda_xxx_finite),
                "M_xxx_norm": finite_probe.m_xxx,
                "ratio_xxx": finite_probe.m_xxx / max(abs(finite_probe.lambda_xxx_finite), 1.0e-30),
                "abs_lambda_xxy_finite": abs(finite_probe.lambda_xxy_finite),
                "M_xxy_norm": finite_probe.m_xxy,
                "ratio_xxy": finite_probe.m_xxy / max(abs(finite_probe.lambda_xxy_finite), 1.0e-30),
                "rigorous_flag": finite_probe.rigorous_flag,
            }
        )

        model_plan: list[tuple[str, SpectralProbe]] = []
        for tail_mode in tail_modes:
            model_plan.append(("norm_bound", probes[tail_mode]))
        model_plan.append(("symmetry_first_order", finite_probe))
        model_plan.append(("signed_xxy_first_order", finite_probe))

        for model, probe in model_plan:
            for shape in cell_shapes:
                radius, values = bisect_cell_size(probe, shape, args.max_radius, args.rxx_samples, model)
                rx, ry = shape_radii(shape, radius)
                rows.append(output_row(model, probe, shape, rx, ry, values))

            for label, rx, ry in target_checks:
                values = jxx_lower(probe, rx, ry, args.rxx_samples, model)
                diagnostics.append(output_row(model, probe, label, rx, ry, values))

    parseval_rows = parseval_decomposition_rows(j_list, modes_d, gx_d, gy_d, modes_n, gx_n, gy_n, weights)
    shell_rows = lambda_xx_shell_contribution_rows(modes_d, gx_d, gy_d, weights, max_j)
    orthogonality_rows = (
        basis_orthogonality_rows(j_list, "D", modes_d, values_d, gx_d, gy_d, weights)
        + basis_orthogonality_rows(j_list, "N", modes_n, values_n, gx_n, gy_n, weights)
    )

    write_outputs(Path(args.outdir), rows, diagnostics, third_rows, parseval_rows, shell_rows, orthogonality_rows)
    print(f"Wrote {args.outdir}/summary.md")
    print(f"Wrote {args.outdir}/max_cell_sizes.csv")
    print(f"Wrote {args.outdir}/diagnostics.csv")
    print(f"Wrote {args.outdir}/third_derivative_diagnostics.csv")
    print(f"Wrote {args.outdir}/parseval_decomposition.csv")
    print(f"Wrote {args.outdir}/lambda_xx_shell_contributions.csv")
    print(f"Wrote {args.outdir}/basis_orthogonality_diagnostics.csv")


if __name__ == "__main__":
    main()
