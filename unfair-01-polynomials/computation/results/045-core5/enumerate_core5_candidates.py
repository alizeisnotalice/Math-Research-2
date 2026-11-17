#!/Users/zzh/.venvs/main/bin/python
"""Exact bounded candidate audit for P=1+a x^2+c x^3+d x^4+x^5.

All filtering is necessary only.  The program writes a compact JSONL first
failure witness for every candidate that passes the discriminant/trace screen.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path


OUTER_POINTS = ((1, 4), (1, 2), (3, 4), (1, 1), (4, 3), (2, 1), (4, 1))


def primitive(coeffs: list[int]) -> list[int]:
    g = 0
    for x in coeffs:
        g = math.gcd(g, abs(x))
    if g > 1:
        coeffs = [x // g for x in coeffs]
    if coeffs[-1] < 0:
        coeffs = [-x for x in coeffs]
    return coeffs


def schur_witness(ascending: list[int]) -> dict | None:
    """Return None iff all roots lie strictly in the unit disk."""
    f = primitive(ascending)
    while len(f) > 1:
        constant, leading = f[0], f[-1]
        if abs(constant) >= abs(leading):
            return {
                "degree": len(f)-1,
                "abs_constant": abs(constant),
                "abs_leading": abs(leading),
            }
        # leading*f - constant*reverse(f), divided by x.  The constant
        # coefficient cancels identically.
        next_f = [leading*f[j] - constant*f[-1-j] for j in range(len(f))]
        assert next_f[0] == 0
        f = primitive(next_f[1:])
    return None


def scaled(coeffs: list[int]) -> list[int]:
    return [h * (2**i) for i, h in enumerate(coeffs)]


def h_coefficients(T: int, N: int, c: int, d: int) -> list[int]:
    return [1, 0, T, 2*c, 2*d+N, 2+T*c, c*c+T*d,
            2*c*d+T, d*d+2*c, 2*d, 1]


def annulus_witness(T: int, N: int, c: int, d: int) -> dict | None:
    h = h_coefficients(T, N, c, d)
    outer = schur_witness(scaled(h))
    if outer is not None:
        return {"orientation": "H(2x)", **outer}
    inner = schur_witness(scaled(list(reversed(h))))
    if inner is not None:
        return {"orientation": "reverse(H)(2x)", **inner}
    return None


def axis_witness(T: int, delta: int, c: int, d: int) -> dict | None:
    for p, q in OUTER_POINTS:
        F = q**5 + c*p**3*q**2 + d*p**4*q + p**5
        K = p**2*q**3
        L = 2*F + T*K
        D = L*L - delta*K*K
        if L <= 0 or D <= 0:
            return {"point": [p, q], "L": L, "product_numerator_x4": D}
    return None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def self_test_schur() -> list[dict]:
    # Root radii respectively 1/2, 2, and 1.  The strict boundary must fail.
    assert schur_witness([-1, 2]) is None       # 2x-1
    assert schur_witness([-2, 1]) is not None   # x-2
    assert schur_witness([-1, 1]) is not None   # x-1
    # Phi_5 has unit-circle roots; after x -> 2x they have radius 1/2.
    assert schur_witness([1, 2, 4, 8, 16]) is None
    return [
        {"polynomial_ascending": [-1, 2], "root_radius": "1/2", "strictly_inside": True},
        {"polynomial_ascending": [-2, 1], "root_radius": "2", "strictly_inside": False},
        {"polynomial_ascending": [-1, 1], "root_radius": "1", "strictly_inside": False},
        {"polynomial_ascending": [1, 2, 4, 8, 16],
         "meaning": "Phi_5(2x), all radii 1/2", "strictly_inside": True},
    ]


def main() -> None:
    schur_tests = self_test_schur()
    root = Path(__file__).resolve().parent
    witness_path = root / "core5_exclusion_witnesses.jsonl"
    delta_path = root / "core5_delta_family_rejections.jsonl"
    survivors: list[list[int]] = []
    counts = {
        "raw_box": 0,
        "delta_trace": 0,
        "positive_axis": 0,
        "annulus_schur": 0,
    }
    axis_by_point = {f"{p}/{q}": 0 for p, q in OUTER_POINTS}
    schur_orientation = {"H(2x)": 0, "reverse(H)(2x)": 0}

    delta_families = 0
    with witness_path.open("w", encoding="utf-8") as witness, \
            delta_path.open("w", encoding="utf-8") as delta_witness:
        for T in range(-19, 20):
            threshold = (20-abs(T))**2
            for N in range(-99, 100):
                # c,d form the innermost raw box.
                counts["raw_box"] += 27*19
                delta = T*T - 4*N
                if delta <= 0 or math.isqrt(delta)**2 == delta or delta >= threshold:
                    delta_families += 1
                    delta_witness.write(json.dumps({
                        "T": T, "N": N, "stage": "delta_trace",
                        "delta": delta, "threshold": threshold,
                        "covered_c_d_tuples": 27*19,
                    }, separators=(",", ":")) + "\n")
                    continue
                for c in range(-13, 14):
                    for d in range(-9, 10):
                        counts["delta_trace"] += 1
                        row = [T, N, c, d]
                        bad_axis = axis_witness(T, delta, c, d)
                        if bad_axis is not None:
                            axis_by_point[f"{bad_axis['point'][0]}/{bad_axis['point'][1]}"] += 1
                            witness.write(json.dumps({
                                "candidate": row, "stage": "positive_axis", "witness": bad_axis
                            }, separators=(",", ":")) + "\n")
                            continue
                        counts["positive_axis"] += 1
                        bad_schur = annulus_witness(T, N, c, d)
                        if bad_schur is not None:
                            schur_orientation[bad_schur["orientation"]] += 1
                            witness.write(json.dumps({
                                "candidate": row, "stage": "annulus_schur", "witness": bad_schur
                            }, separators=(",", ":")) + "\n")
                            continue
                        counts["annulus_schur"] += 1
                        survivors.append(row)

    final_path = root / "core5_survivors.json"
    final_path.write_text(json.dumps(survivors, indent=2) + "\n", encoding="utf-8")
    result = {
        "parameters": {
            "T": [-19, 19], "N": [-99, 99], "c": [-13, 13], "d": [-9, 9],
            "positive_axis_points": [list(x) for x in OUTER_POINTS],
            "delta_conditions": [
                "delta=T^2-4N>0",
                "delta is nonsquare",
                "delta<(20-|T|)^2",
            ],
        },
        "counts": counts,
        "axis_first_failure_counts": axis_by_point,
        "schur_first_failure_counts": schur_orientation,
        "survivor_count": len(survivors),
        "survivors_sha256": sha256(final_path),
        "witnesses_sha256": sha256(witness_path),
        "delta_family_witnesses_sha256": sha256(delta_path),
        "coverage": {
            "delta_rejected_families": delta_families,
            "delta_rejected_tuples": delta_families * 27 * 19,
            "post_delta_tuples": counts["delta_trace"],
            "post_delta_partition_check": (
                (counts["delta_trace"] - counts["positive_axis"])
                + (counts["positive_axis"] - counts["annulus_schur"])
                + counts["annulus_schur"]
                == counts["delta_trace"]
            ),
        },
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
        },
        "schur_semantics_self_tests": schur_tests,
        "scope": (
            "Every screen is a necessary condition. Surviving tuples are not "
            "claimed to define P, Pbar, a complete-cross model, or a counterexample."
        ),
    }
    (root / "core5_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
