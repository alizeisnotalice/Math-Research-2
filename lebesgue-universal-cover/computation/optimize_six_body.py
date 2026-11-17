"""Floating local search for hull(D,B3,B5,B7,B9,P11).

This is deliberately a numerical probe: Reuleaux boundaries are sampled, so
the objective is an inscribed-polygon hull area, not a certified six-body
value.  The parameterization fixes the disk and B3 orientation using global
rigid-motion freedom; it still has 14 relative degrees of freedom.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.spatial import ConvexHull


OUT = Path(__file__).with_name("six_body_numerical_result.json")
SEED = 20261117
DEADLINE_SECONDS = 155.0
ARC_SAMPLES = 48
DISK_SAMPLES = 384

# Tomar's opt11 placement: B3 translation; then (angle, translation) for
# B5, B7, B9.  The new P11 uses its best floating relaxation placement from
# round 5 as the first start.
X_FIVE = np.array([
    -0.0021161707485670258, -0.002171588117597346,
    0.03368534128826296, 0.025084694928075076, -0.0038913081911800527,
    0.708816831902889, 0.0035625741962694224, -0.012629515020220292,
    0.3490932487151976, -0.0006351709996987754, -0.0006510120748306202,
])
P11_START = np.array([0.14914631789769725, 0.0018446783402523545, 0.0021145725048262935])
NS = (3, 5, 7, 9)


class TimeLimit(Exception):
    pass


def u(a: np.ndarray) -> np.ndarray:
    return np.column_stack((np.cos(a), np.sin(a)))


def reuleaux_points(n: int, rho: float, tx: float, ty: float) -> np.ndarray:
    """Boundary samples from the exact radius-one circular arcs."""
    R = 1.0 / (2.0 * math.cos(math.pi / (2.0 * n)))
    out = []
    for ell in range(n):
        phi = rho + 2.0 * math.pi * ell / n
        center = np.array((tx + R * math.cos(phi), ty + R * math.sin(phi)))
        a = np.linspace(
            phi + math.pi - math.pi / (2.0 * n),
            phi + math.pi + math.pi / (2.0 * n),
            ARC_SAMPLES + 1,
        )
        out.append(center + u(a))
    return np.vstack(out)


def all_points(x: np.ndarray) -> np.ndarray:
    assert x.shape == (14,)
    disk = 0.5 * u(np.linspace(0.0, 2.0 * math.pi, DISK_SAMPLES, endpoint=False))
    groups = [disk, reuleaux_points(3, 0.0, x[0], x[1])]
    k = 2
    for n in (5, 7, 9):
        groups.append(reuleaux_points(n, x[k], x[k + 1], x[k + 2]))
        k += 3
    assert k == 11
    R11 = 1.0 / (2.0 * math.cos(math.pi / 22.0))
    a11 = x[11] + 2.0 * math.pi * np.arange(11) / 11.0
    groups.append(np.array((x[12], x[13])) + R11 * u(a11))
    return np.vstack(groups)


def hull_area(x: np.ndarray) -> float:
    return float(ConvexHull(all_points(x)).volume)


def main() -> None:
    started = time.monotonic()
    calls = 0
    best_x = np.r_[X_FIVE, P11_START]
    best_area = hull_area(best_x)

    bounds = [(-0.1, 0.1), (-0.1, 0.1)]
    for n in (5, 7, 9):
        bounds.extend([(0.0, 2.0 * math.pi / n), (-0.1, 0.1), (-0.1, 0.1)])
    bounds.extend([(0.0, 2.0 * math.pi / 11.0), (-0.1, 0.1), (-0.1, 0.1)])
    assert len(bounds) == 14

    def objective(x: np.ndarray) -> float:
        nonlocal calls, best_x, best_area
        calls += 1
        if time.monotonic() - started > DEADLINE_SECONDS:
            raise TimeLimit
        val = hull_area(x)
        if val < best_area:
            best_area, best_x = val, np.array(x, copy=True)
        return val

    phases = []
    # P11 alone is optimized first, against the five-body source placement.
    try:
        sol = minimize(
            lambda z: objective(np.r_[X_FIVE, z]), P11_START,
            method="Powell", bounds=bounds[11:],
            options={"maxiter": 450, "xtol": 2e-9, "ftol": 2e-11},
        )
        phases.append({"phase": "P11-only", "success": bool(sol.success), "message": str(sol.message), "fun": float(sol.fun), "nit": int(sol.nit), "nfev": int(sol.nfev)})
        x0 = np.r_[X_FIVE, sol.x]
    except TimeLimit:
        phases.append({"phase": "P11-only", "stopped": "deadline"})
        x0 = best_x

    # Two local joint starts: the first continuation and one tiny deterministic
    # perturbation.  They do not claim global search coverage.
    rng = np.random.default_rng(SEED)
    starts = [x0, np.clip(x0 + rng.normal(0.0, 0.002, size=14), [b[0] for b in bounds], [b[1] for b in bounds])]
    for i, start in enumerate(starts):
        try:
            sol = minimize(
                objective, start, method="Powell", bounds=bounds,
                options={"maxiter": 1100, "xtol": 2e-8, "ftol": 2e-10},
            )
            phases.append({"phase": f"joint-{i}", "success": bool(sol.success), "message": str(sol.message), "fun": float(sol.fun), "nit": int(sol.nit), "nfev": int(sol.nfev)})
        except TimeLimit:
            phases.append({"phase": f"joint-{i}", "stopped": "deadline"})
            break

    # Finer sampling independently evaluates the best returned placement.
    global ARC_SAMPLES, DISK_SAMPLES
    coarse_area = hull_area(best_x)
    ARC_SAMPLES, DISK_SAMPLES = 192, 1536
    fine_area = hull_area(best_x)
    elapsed = time.monotonic() - started
    payload = {
        "status": "completed",
        "scope": "floating local search only; sampled boundary hull and Powell do not certify a six-body minimum or a lower bound",
        "seed": SEED,
        "elapsed_seconds": elapsed,
        "objective_calls": calls,
        "parameterization": "disk fixed at origin; B3 angle fixed to zero; x=[t3x,t3y,rho5,t5x,t5y,rho7,t7x,t7y,rho9,t9x,t9y,rho11,t11x,t11y]",
        "degrees_of_freedom": 14,
        "source_five_body_params": X_FIVE.tolist(),
        "best_parameters": best_x.tolist(),
        "coarse_sample_area": coarse_area,
        "fine_sample_area": fine_area,
        "coarse_sampling": {"arc_segments_per_Reuleaux_arc": 48, "disk_vertices": 384},
        "fine_sampling": {"arc_segments_per_Reuleaux_arc": 192, "disk_vertices": 1536},
        "phases": phases,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
