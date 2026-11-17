"""Floating screen of constant-width Fourier bodies against a certified cover.

For h=1/2+a cos(3t)+b cos(5t+phi), h+h'' >= 0 follows from
8|a|+24|b| <= 1/2.  This script tests only 64 sampled boundary points and a
finite angular grid in an outer half-plane relaxation, so a positive value is
only a candidate for later interval certification.
"""
from __future__ import annotations
import hashlib, json, math, platform, sys, time
from fractions import Fraction
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linprog

WORK = Path("/Volumes/data/CodexProjects/T2/.local/workspaces/lebesgue-main/lebesgue-universal-cover/computation")
P10, P7 = WORK / "010-seven-body-certificate.json", WORK / "007-six-body-certificate.json"
CERT = P10 if P10.is_file() and P10.stat().st_size else P7
OUT = Path(__file__).with_name("fourier_width_probe_result.json")
SEED = 20261117
BOUNDARY_POINTS = 64
COARSE_ANGLES = 72
FINE_ANGLES = 360
DEADLINE_SECONDS = 150.0


class Stop(Exception): pass


def body(a: float, b: float, phi: float, reflected: bool) -> np.ndarray:
    t = 2.0 * math.pi * np.arange(BOUNDARY_POINTS) / BOUNDARY_POINTS
    h = 0.5 + a * np.cos(3.0 * t) + b * np.cos(5.0 * t + phi)
    hp = -3.0 * a * np.sin(3.0 * t) - 5.0 * b * np.sin(5.0 * t + phi)
    pts = np.column_stack((h * np.cos(t) - hp * np.sin(t), h * np.sin(t) + hp * np.cos(t)))
    if reflected: pts[:, 1] *= -1.0
    return pts


def main() -> None:
    began = time.monotonic()
    cert_bytes = CERT.read_bytes()
    source = json.loads(cert_bytes)
    V = np.array(source["hull_grid_points"], dtype=float) / source["coordinate_denominator"]
    edge = np.roll(V, -1, axis=0) - V
    N = np.column_stack((edge[:, 1], -edge[:, 0]))
    N /= np.linalg.norm(N, axis=1)[:, None]
    eps = float(Fraction(source["epsilon"]))
    rhs_cover = (N * V).sum(axis=1) + eps
    A = np.column_stack((N, -np.ones(len(N))))
    calls = 0

    def min_slack(a: float, b: float, phi: float, n_angles: int) -> dict:
        nonlocal calls
        best = None
        for branch in (False, True):
            P = body(a, b, phi, branch)
            for j in range(n_angles):
                if time.monotonic() - began > DEADLINE_SECONDS:
                    raise Stop
                theta = 2.0 * math.pi * j / n_angles
                c, s = math.cos(theta), math.sin(theta)
                Q = P @ np.array(((c, s), (-s, c)))
                h = (N @ Q.T).max(axis=1)
                sol = linprog([0.0, 0.0, 1.0], A_ub=A, b_ub=rhs_cover - h, bounds=[(None, None)] * 3, method="highs")
                calls += 1
                if not sol.success:
                    raise RuntimeError(sol.message)
                row = {"slack": float(sol.fun), "theta": theta, "translation": sol.x[:2].tolist(), "reflected": branch}
                if best is None or row["slack"] < best["slack"]:
                    best = row
        return best

    # q=16a and r=48b make 8a+24b=(q+r)/2.  Every row below has q+r=.96,
    # giving a positive 0.02 curvature-radius reserve before floating error.
    qr = [(0.96, 0.0), (0.84, 0.12), (0.70, 0.26), (0.54, 0.42), (0.38, 0.58), (0.20, 0.76), (0.0, 0.96)]
    candidates = []
    for q, r in qr:
        for k in range(12):
            candidates.append((q / 16.0, r / 48.0, 2.0 * math.pi * k / 12.0, q, r))
    coarse = []
    try:
        for a, b, phi, q, r in candidates:
            coarse.append({"a": a, "b": b, "phi": phi, "q": q, "r": r, "curvature_margin": 0.5 - 8.0 * abs(a) - 24.0 * abs(b), "coarse": min_slack(a, b, phi, COARSE_ANGLES)})
    except Stop:
        pass
    coarse.sort(key=lambda row: row["coarse"]["slack"], reverse=True)
    refined = []
    try:
        for row in coarse[:5]:
            refined.append({**{k: row[k] for k in ("a", "b", "phi", "q", "r", "curvature_margin")}, "fine": min_slack(row["a"], row["b"], row["phi"], FINE_ANGLES)})
    except Stop:
        pass
    refined.sort(key=lambda row: row["fine"]["slack"], reverse=True)
    samples = []
    if refined:
        best = refined[0]
        samples = body(best["a"], best["b"], best["phi"], best["fine"]["reflected"]).tolist()
    result = {
        "status": "completed" if len(coarse) == len(candidates) else "deadline_before_full_grid",
        "scope": "floating finite-angle LP over an outer half-plane relaxation and 64 boundary samples; positive slack is only a later-certification candidate",
        "cover_certificate": str(CERT),
        "cover_certificate_sha256": hashlib.sha256(cert_bytes).hexdigest(),
        "cover_vertices": int(len(V)),
        "cover_epsilon": source["epsilon"],
        "formula": "h(t)=1/2+a cos(3t)+b cos(5t+phi); x=h*u+hprime*u_perp",
        "constraint": "8|a|+24|b|<=1/2",
        "boundary_points": BOUNDARY_POINTS,
        "coarse_orientation_grid": COARSE_ANGLES,
        "fine_orientation_grid": FINE_ANGLES,
        "branches": ["identity", "reflection in x-axis before rotation"],
        "candidate_family": "q=16a,r=48b with q+r=0.96; phi=2pi*k/12; nonnegative a,b subfamily only",
        "coarse_results": coarse,
        "refined_results": refined,
        "best_boundary_samples": samples,
        "lp_calls": calls,
        "elapsed_seconds": time.monotonic() - began,
        "seed": SEED,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__": main()
