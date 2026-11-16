#!/usr/bin/env python3
"""Minimal Fraction kernel for the R81 analytic shared-side payment."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def plain(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--producer-family", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    start = time.perf_counter()
    cfg = json.loads(args.input.read_text())
    assert sha(args.producer_family) == cfg["producer_family_sha256"]
    raw = json.loads(args.producer_family.read_text())
    w, e, m, d, D = map(Q, [cfg[k] for k in
                          ["w", "e", "old_corner_mass", "near_distance", "far_distance"]])
    tlo, thi = map(Q, cfg["parameter_interval"])
    assert w == 1-e and m == Q(1, 4)
    assert raw["S_indices"] == [0, 1, 2, 3] and raw["threshold"] == "1"
    assert list(map(Q, raw["requested_t_interval"])) == [tlo, thi]

    def turn(p):
        return (w-p[1], p[0])

    positions = [(Q(0), Q(0)), (w, Q(0)), (w, w), (Q(0), w)]
    coeffs = [(Q(0), m)]*4
    q, z = (-d, Q(0)), (Q(0), -D)
    for i in range(4):
        positions.extend([q, z])
        coeffs.extend([(Q(0), m), (Q(1), Q(0))] if i % 2 == 0
                      else [(Q(-1), Q(1, 2)), (Q(0), m)])
        q, z = turn(q), turn(z)
    assert positions == [tuple(map(Q, p)) for p in raw["positions"]]
    assert coeffs == [tuple(map(Q, p)) for p in raw["mass_affine_slope_intercept"]]
    assert len(set(positions)) == 12
    side_affine = []
    for i in range(4):
        near = coeffs[4+2*i]
        far = coeffs[5+2*((i-1) % 4)]
        combined = (near[0]+far[0], near[1]+far[1])
        assert combined == (0, Q(1, 2))
        side_affine.append(combined)
    assert sum(s for s, intercept in coeffs[4:]) == 0
    assert sum(intercept for s, intercept in coeffs[4:]) == 2
    assert Q(1, 2)-thi == Q(1, 1000) > 0

    # Directed bound on the even near component; square comparison is exact.
    a_lower = Q(176776, 500000)
    a_upper = Q(176777, 500000)
    assert a_lower*a_lower < Q(1, 8) < a_upper*a_upper
    assert e < d-a_upper < m < a_lower < a_upper < d
    strong_upper = Q(182323, 500000)-Q(171, 200)*a_lower
    assert strong_upper < Q(63, 1000)
    # Odd-corner G(a): two analytic cases, not sampled t values.
    C0 = m/2+d*m-m*e+e*e
    assert C0 == Q(213501, 1000000)
    inactive_mass_bound = 4*(d-m)**2
    inactive_upper = C0+inactive_mass_bound/2-(d+m+e)*m
    assert inactive_mass_bound == Q(441, 10000)
    assert inactive_upper == Q(84051, 1000000) < Q(1, 10)
    A, B = d+m, d+e
    C = C0-d*e
    assert (A, B, C) == (Q(605, 1000), Q(356, 1000), Q(213146, 1000000))
    reverse_cauchy_gap = A*A-B*B-Q(48, 100)**2
    assert reverse_cauchy_gap > 0
    active_upper = C-m*Q(48, 100)
    assert active_upper == Q(93146, 1000000) < Q(1, 10)

    # Joint activity threshold, derivative cap and complete hinge branches.
    B0 = D*D-m
    kappa = Q(1, 2)-(D+2*e)/4
    delta = Q(1, 2)-B0
    etamax = thi-tlo
    assert B0 == Q(19241, 40000)
    assert kappa == Q(1143, 4000)
    assert delta == Q(759, 40000) and 0 < delta < etamax
    # On eta<=delta the sum is 4delta; on eta>=delta it is 2(delta+eta).
    branch1 = 4*delta
    branch2max = 2*(delta+etamax)
    assert branch1 <= branch2max == Q(10719, 20000)
    joint_upper = kappa*branch2max
    assert joint_upper < Q(31, 200)
    near_upper = 2*Q(63, 1000)+2*Q(1, 10)
    corner_upper = near_upper+Q(31, 200)
    noncorner_area = 3-3*e*e
    complete_upper = noncorner_area+corner_upper
    quota = Q(1, 2)+3*e*e
    assert corner_upper == Q(481, 1000) < quota
    assert complete_upper == Q(3480997, 1000000) < Q(7, 2)
    # Source orderings hold throughout the allowed ranges.
    assert Q(499, 4000) < (D-Q(1, 2))**2
    assert Q(749, 4000) < Q(433, 1000)**2
    assert D-Q(433, 1000) > a_upper
    assert D-Q(1, 2) > a_upper

    # Uncapped m=1 affine-cost obstruction: exact radical comparison only.
    upper101 = Q(443, 44)
    assert upper101**2 > 101
    affine_threshold = Q(1, 4)+Q(1, 100)/8
    assert affine_threshold == Q(201, 800)
    assert (Q(161)-11*upper101)/200 == affine_threshold
    assert Q(101, 100) > 1  # Explicit failure of global unit cap at old mass 1.

    result = {"status": "PASS", "scope": cfg["scope"],
              "family_source_binding": {"positions": positions, "mass_affine": coeffs,
                                        "four_shared_side_mass_affine": side_affine},
              "analytic_constants": {"even_near_upper": strong_upper,
                "odd_inactive_upper": inactive_upper, "odd_active_upper": active_upper,
                "reverse_Cauchy_squared_margin": reverse_cauchy_gap,
                "joint_threshold_B0": B0, "joint_derivative_upper": kappa,
                "hinge_constant_branch": branch1, "hinge_global_upper": branch2max,
                "joint_sum_upper": joint_upper, "near_sum_strict_upper": near_upper,
                "corner_sum_strict_upper": corner_upper, "noncorner_container_area": noncorner_area,
                "complete_deletion_strict_upper": complete_upper,
                "target_margin": Q(7, 2)-complete_upper},
              "uncapped_affine_counterexample": {
                "exact_area": "(161-11sqrt(101))/200", "strictly_above": affine_threshold,
                "sqrt101_strict_upper": upper101, "global_unit_cap_fails": True},
              "receipt": {"pid": os.getpid(), "python": sys.version,
                "platform": platform.platform(), "working_directory": str(Path.cwd()),
                "actual_command": [sys.executable]+sys.argv, "input_sha256": sha(args.input),
                "script_sha256": sha(__file__), "producer_family_sha256": sha(args.producer_family),
                "seconds": time.perf_counter()-start, "random_seed": None,
                "tolerance": "exact Fraction; none", "area_or_parameter_search_performed": False,
                "cap_polytope_rerun": False}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plain(result), indent=2)+"\n")
    print(json.dumps({"status": "PASS", "corner_sum_strict_upper": str(corner_upper),
                      "complete_deletion_strict_upper": str(complete_upper),
                      "joint_sum_upper": str(joint_upper)}))


if __name__ == "__main__":
    main()
