#!/usr/bin/env python3
"""Minimal exact scalar replay for the R79 count/envelope proof; no area code."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameters(n):
    return Q(1, n*n), Q(n+1, n*n), Q(n*n, (n+1)**2)


def ceilq(a):
    return -((-a.numerator)//a.denominator)


def soft(n, X, Y, r):
    e, h, c = parameters(n)
    return c*(r-X+h)*(r-Y+h)


def f(n, X, Y, r):
    return 4*r*r-soft(n, X, Y, r)


def minimum(n, X, Y, R):
    """Open-radius existential test uses B>phi, after R<1/2 gate."""
    if R >= Q(1, 2):
        return None
    e, h, c = parameters(n)
    r0 = c*(2*h-X-Y)/(2*(4-c))
    rs = min(Q(1, 2), max(R, r0))
    return rs, f(n, X, Y, rs)


def plain(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    return value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--actual-input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    started = time.perf_counter()
    data = json.loads(args.input.read_text())
    assert sha(args.actual_input) == data["actual_input_sha256"]
    raw = json.loads(args.actual_input.read_text())
    positions = [tuple(map(Q, p)) for p in raw["positions"]]
    masses = list(map(Q, raw["masses"]))
    old = [(Q(21*i, 400), Q(21*j, 400)) for i in range(20) for j in range(20)]
    assert positions[:400] == old
    assert positions[400:] == [(Q(-71, 200), Q(0)), (Q(0), Q(-171, 200))]
    assert masses == [Q(1, 400)]*400+[Q(1, 4), Q(1, 2)]
    assert raw["S_indices"] == list(range(400))

    rows = []
    integer_event_count = 0
    for n in data["n_values"]:
        e, h, c = parameters(n)
        epsilon = Q(n*n+n-1, n*n*(n+1))
        assert epsilon == (1-2*e)/(n+1)+e
        # Endpoints and integer source events are deliberate, not random sampling.
        for X in [e, e+h/4]:
            for Y in [e, e+h/3]:
                for r in [max(X, Y)+h/2, X+h, Y+h, Q(49, 100)]:
                    if not max(X, Y) < r < Q(1, 2):
                        continue
                    a, b = (r-X)/h, (r-Y)/h
                    nx = sum(Q(i)*h < r-X for i in range(n))
                    ny = sum(Q(j)*h < r-Y for j in range(n))
                    assert (nx, ny) == (ceilq(a), ceilq(b))
                    integer_event_count += int(a.denominator == 1 or b.denominator == 1)
                    actual = e*nx*ny
                    lo, up = c*(r-X)*(r-Y), soft(n, X, Y, r)
                    local_error = c*h*(2*r-X-Y+h)
                    tx, ty = Q(nx)-a, Q(ny)-b
                    assert up-actual == e*(a*(1-ty)+b*(1-tx)+1-tx*ty)
                    assert lo <= actual < up
                    assert 0 < up-actual <= local_error < epsilon
                    expanded = (4-c)*r*r+c*(X+Y-2*h)*r-c*(X-h)*(Y-h)
                    assert f(n, X, Y, r) == expanded
                    base = (r-X)*(r-Y)
                    assert 0 <= up-base <= (1-c)/4+Q(1, n+1)+e
                    rows.append({"n": n, "X": X, "Y": Y, "r": r,
                                 "actual": actual, "soft": up, "error_bound": local_error})
    assert integer_event_count > 0

    # Real-domain convex branches, strict equality, B=0, and invalid interval.
    branches = []
    for n, X, Y, R in [(2, Q(1, 4), Q(1, 4), Q(1, 4)),
                        (5, Q(1, 25), Q(1, 25), Q(1, 25)),
                        (20, Q(2, 5), Q(39, 100), Q(93, 200))]:
        e, h, c = parameters(n)
        rs, phi = minimum(n, X, Y, R)
        r0 = c*(2*h-X-Y)/(2*(4-c))
        assert rs == max(R, r0)
        assert r0 <= Q(n, 3*n*n+8*n+4) < Q(1, 2)
        assert not phi > phi  # Strict test at B=phi must reject.
        assert phi+Q(1, 1000000) > phi
        for r in [R, Q(1, 2), rs, (R+Q(1, 2))/2]:
            assert f(n, X, Y, r) >= phi
        branches.append({"n": n, "radius_minimizer": rs, "minimum": phi,
                         "B_zero_accepted": Q(0) > phi,
                         "vertex_inside": R < r0 < Q(1, 2)})
    assert any(t["vertex_inside"] for t in branches)
    assert any(t["B_zero_accepted"] for t in branches)
    assert any(not t["B_zero_accepted"] for t in branches)
    assert minimum(20, Q(2, 5), Q(39, 100), Q(1, 2)) is None
    assert minimum(20, Q(2, 5), Q(39, 100), Q(3, 5)) is None

    n = 20
    X, Y = map(Q, data["comparison_point"])
    qs, bs = positions[400:], masses[400:]
    comparisons = []
    for mask in range(4):
        selected = [i for i in range(2) if mask & (1 << i)]
        R = max([X, Y]+[max(abs(qs[i][0]+X), abs(qs[i][1]+Y)) for i in selected])
        B = sum((bs[i] for i in selected), Q(0))
        rs, phi = minimum(n, X, Y, R)
        assert rs == R and phi > B
        comparisons.append({"mask": mask, "R": R, "mass": B, "minimum": phi,
                            "rejection_margin": phi-B})
    assert comparisons[0]["minimum"] == Q(16, 25)-Q(1, 336)
    assert comparisons[3]["minimum"] == Q(8649, 10000)-Q(799, 58800)
    residual = []
    for mask in [1, 2, 3]:
        selected = [i for i in range(2) if mask & (1 << i)]
        d = max(max(abs(qs[i][0]+X), abs(qs[i][1]+Y)) for i in selected)
        B = sum((bs[i] for i in selected), Q(0))
        assert 4*d*d > B
        residual.append({"mask": mask, "distance": d, "mass": B,
                         "residual_rejection_margin": 4*d*d-B})
    assert comparisons[3]["R"] < Q(1, 2)
    assert Q(1, 4)+sum(bs) == 1  # Total-local-mass proxy radius is 1/2.
    special = data["raw_mass_exceeds_proxy"]
    excess = soft(special["n"], Q(special["X"]), Q(special["Y"]), Q(special["r"]))-Q(1, 4)
    assert excess > 0

    result = {"status": "PASS", "scope": data["scope"],
              "open_count_fixtures": rows, "open_count_fixture_count": len(rows),
              "integer_event_fixture_count": integer_event_count,
              "convex_branches": branches, "comparison_point": [X, Y],
              "all_four_subset_rejections": comparisons, "complete_residual_rejections": residual,
              "total_mass_proxy_accepts": True, "raw_soft_exceeds_total_mass_by": excess,
              "receipt": {"python": sys.version, "platform": platform.platform(),
                          "actual_command": [sys.executable]+sys.argv,
                          "working_directory": str(Path.cwd()),
                          "input_sha256": sha(args.input), "script_sha256": sha(__file__),
                          "actual_input_sha256": sha(args.actual_input),
                          "wall_seconds": time.perf_counter()-started,
                          "exit_status": "successful completion; command exit code recorded by caller",
                          "random_seed": None, "tolerance": "exact rational; none"}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plain(result), indent=2)+"\n")
    print(json.dumps({"status": "PASS", "count_fixtures": len(rows),
                      "integer_events": integer_event_count,
                      "minimum_rejection_margin": str(min(x["rejection_margin"] for x in comparisons))}))


if __name__ == "__main__":
    main()
