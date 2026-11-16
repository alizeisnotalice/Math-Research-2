#!/usr/bin/env python3
"""R82 exact symbolic Laurent/Fraction kernel, without numerical sampling."""
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
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    return x


def add(*ps):
    out = {}
    for p in ps:
        for power, coefficient in p.items():
            out[power] = out.get(power, Q(0))+coefficient
    return {power: coefficient for power, coefficient in out.items() if coefficient}


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i+j] = out.get(i+j, Q(0))+a*b
    return {power: coefficient for power, coefficient in out.items() if coefficient}


def const(q):
    return {0: Q(q)} if q else {}


def scale(q, p):
    return mul(const(q), p)


def sub(p, q):
    return add(p, scale(-1, q))


def derivative(p):
    return {power-1: coefficient*power for power, coefficient in p.items() if power}


def value(p, s):
    return sum((coefficient*s**power for power, coefficient in p.items()), Q(0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    start = time.perf_counter()
    cfg = json.loads(args.input.read_text())
    m, v, d, e = map(Q, [cfg[k] for k in ["old_mass", "old_mass_radius", "near_distance", "crop_inset"]])
    amin, amax = map(Q, cfg["outside_mass_interval"])
    target, comparison, s0 = map(Q, [cfg[k] for k in ["target", "intermediate_comparison", "high_branch_anchor_s"]])
    assert m == 4*v*v == Q(1, 4) and amin == 0 and amax == Q(1, 2)
    S = {1: Q(1)}
    U = {1: Q(1, 2), -1: v*v/2}
    Z = {1: Q(1, 2), -1: -v*v/2}
    assert sub(mul(U, U), mul(Z, Z)) == const(v*v)
    assert add(U, Z) == S
    C0, L = m/2+d*v-v*e+e*e, d+v+e
    full_before_cover = add(const(C0), scale(2, mul(Z, Z)), scale(-L, U))
    contact_term = mul(sub(S, const(d)), sub(Z, const(e)))
    G2 = sub(full_before_cover, contact_term)
    G3 = sub(mul(sub(U, const(e)), sub(U, const(e))), contact_term)
    assert G2 == {0: Q(22737, 125000), 1: -Q(249, 2000), -1: -Q(961, 32000), -2: Q(1, 512)}
    assert G3 == {2: -Q(1, 4), 1: Q(71, 400), 0: Q(31073, 500000),
                  -1: -Q(357, 32000), -2: Q(1, 1024)}
    # Both full-area formulas match when d-u-e=0, identically before imposing it.
    assert sub(G3, G2) == mul(sub(const(d-e), U), sub(U, const(v)))

    zero_cost = (v-e)**2
    assert value(U, v) == v and value(Z, v) == 0
    assert value(full_before_cover, v) == zero_cost == Q(62001, 1000000) < target
    assert v+2*e < d  # z<=e implies s<=v+2z<=v+2e<d.
    assert v > d/2
    umax_up, zmax_up = map(Q, [cfg["sqrt_u_max_upper"], cfg["sqrt_z_max_upper"]])
    assert umax_up**2 > (m+amax)/4
    assert zmax_up**2 > amax/4
    assert umax_up+zmax_up < 1
    assert zmax_up < d-e and umax_up < Q(1, 2)

    # Whole low-contact branch: monotonic derivative and exact rational endpoint.
    monotone_margin = 4-L/v
    assert monotone_margin == Q(197, 125) > 0
    u_contact, z_contact = value(U, d), value(Z, d)
    assert u_contact < d-e and z_contact > e
    contact_value = value(full_before_cover, d)
    assert contact_value == value(G2, d) == Q(691638777, 10082000000) < target

    # Whole intermediate branch: exact polynomial Taylor positivity.
    P = mul({2: Q(1)}, sub(const(comparison), G2))
    assert P == {3: Q(249, 2000), 2: -Q(13987, 125000), 1: Q(961, 32000), 0: -Q(1, 512)}
    P0, P1, P2 = value(P, d), value(derivative(P), d), value(derivative(derivative(P)), d)
    assert (P0, P1, P2) == (Q(14101223, 80000000000), -Q(937829, 400000000), Q(41393, 1000000))
    assert value(derivative(derivative(derivative(P))), d)/6 == P[3] > 0
    K = P0-P1*P1/(2*P2)
    assert K == Q(291048892263, 2649152000000000) > comparison-target == Q(1, 10000)
    # Thus P>=K on s>=d, and P/s²>K when s<1. No s grid is used.

    # Whole high-cover branch: anchor lies before coverage, derivative is negative.
    assert value(U, s0) == Q(169, 480) < d-e == Q(177, 500)
    expected_derivative = {-1: Q(0), 1: -Q(1, 2), 0: d/2,
                           -2: (d+2*e)*v*v/2, -3: -v**4/2}
    expected_derivative = {k: x for k, x in expected_derivative.items() if x}
    assert derivative(G3) == expected_derivative
    high_derivative_upper = (d-s0)/2+(d+2*e)*v*v/(2*s0*s0)
    assert high_derivative_upper == -Q(1757, 19200) < 0
    high_value_upper = value(G3, s0)
    assert high_value_upper == Q(9038149, 144000000) < target

    result = {"status": "PASS", "scope": cfg["scope"],
              "exact_symbolic_laurent_coefficients": {"u": U, "z": Z,
                "full_before_cover": full_before_cover, "G2": G2, "G3": G3, "P": P},
              "zero_and_low_branch": {"zero_mass_exact_cost": zero_cost,
                "z_below_crop_implies_no_contact": v+2*e < d,
                "monotonic_derivative_margin": monotone_margin,
                "contact_u": u_contact, "contact_z": z_contact, "contact_cost": contact_value},
              "intermediate_branch": {"P_at_d": P0, "Pprime_at_d": P1,
                "Psecond_at_d": P2, "positive_cubic_remainder_coefficient": P[3],
                "Taylor_global_lower": K, "comparison_margin": comparison-target,
                "s_less_than_one_from_root_bounds": True},
              "high_branch": {"anchor_s": s0, "anchor_u": value(U, s0),
                "coverage_u_threshold": d-e, "derivative_uniform_upper": high_derivative_upper,
                "cost_at_anchor": high_value_upper},
              "uniform_actual_near_area_strict_upper": target,
              "receipt": {"pid": os.getpid(), "python": sys.version,
                "platform": platform.platform(), "actual_command": [sys.executable]+sys.argv,
                "working_directory": str(Path.cwd()), "input_sha256": sha(args.input),
                "script_sha256": sha(__file__), "seconds": time.perf_counter()-start,
                "random_seed": None, "tolerance": "exact Fraction; none",
                "mass_sample_count": 0, "parameter_subdivision_performed": False,
                "area_union_integration_performed": False}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plain(result), indent=2)+"\n")
    print(json.dumps({"status": "PASS", "uniform_actual_near_area_strict_upper": str(target),
                      "Taylor_global_lower": str(K), "mass_samples": 0}))


if __name__ == "__main__":
    main()
