"""Select numerical degree-56 root-allocation targets using exact prime tests."""

import argparse
from collections import Counter
import json
from math import lcm
from pathlib import Path

import numpy as np
import sympy as sp


X = sp.Symbol("x")
H = (1, -4, 6, -7, 4, -2, 1)
PRIMES = (2, 3, 7, 11, 163, 1279)
DEGREE = 56
FREQUENCIES = 64
RADIUS = 1.031


def factor_period(coefficients, prime):
    f = sp.Poly.from_list(coefficients, X, modulus=prime)
    _, parts = f.factor_list()
    degrees = sorted((g.degree(), exponent) for g, exponent in parts)
    assert sum(d*e for d, e in degrees) == f.degree()
    return lcm(*(d for d, _ in degrees)), degrees


def root_packet(c, points):
    roots = np.polynomial.polynomial.polyroots(c)
    reconstruction = np.polynomial.polynomial.polyfromroots(roots)
    error = float(np.max(abs(reconstruction - c)))
    if error > 1e-7:
        return None, "root_reconstruction_tolerance"
    unused = set(range(DEGREE))
    groups = []
    while unused:
        i = min(unused)
        unused.remove(i)
        z = roots[i]
        if abs(z.imag) < 1e-8:
            groups.append([complex(z.real)])
        else:
            if not unused:
                return None, "unpaired_root"
            j = min(unused, key=lambda j: abs(roots[j]-z.conjugate()))
            if abs(roots[j]-z.conjugate()) > 1e-8:
                return None, "unpaired_root"
            groups.append([z, roots[j]])
            unused.remove(j)
    if not 1 < len(groups) < 32:
        return None, "native_group_capacity"
    logs = np.array([np.log(np.prod([(points-z)/(-z) for z in group], axis=0)) for group in groups])
    if not np.isfinite(logs).all():
        return None, "nonfinite_log"
    return {
        "C": c.tolist(), "roots": [[[float(z.real), float(z.imag)] for z in group] for group in groups],
        "LR": logs.real.tolist(), "LI": logs.imag.tolist(), "groups": len(groups),
        "masks": 1 << (len(groups)-1), "root_product_error": error,
    }, None


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tasks", type=int, default=12)
    parser.add_argument("--draw-cap", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261118056)
    args = parser.parse_args()
    assert Path("/Volumes/data").is_mount()
    periods = {p: factor_period(H, p)[0] for p in PRIMES}
    assert list(periods.values()) == [6, 4, 6, 5, 2, 3]
    rng = np.random.default_rng(args.seed)
    points = RADIUS * np.exp(2j*np.pi*np.arange(FREQUENCIES)/FREQUENCIES)
    tasks, rejections, visited = [], [], set()
    for draw in range(args.draw_cap):
        c = rng.integers(0, 2, DEGREE+1)
        c[:12], c[-1] = [1]*7+[0]+[1]*4, 1
        key = tuple(map(int, c))
        if key in visited:
            continue
        visited.add(key)
        local_periods = {}
        for prime, needed in periods.items():
            period, degrees = factor_period(key[::-1], prime)
            if period % needed:
                rejections.append({"draw": draw, "C": key, "kind": "exact_residue_period",
                                   "prime": prime, "required": needed, "degrees": degrees})
                break
            local_periods[prime] = {"period": period, "degrees": degrees}
        else:
            data, failure = root_packet(c, points)
            if failure:
                rejections.append({"draw": draw, "C": key, "kind": failure})
                continue
            tasks.append({"id": len(tasks), "draw": draw, "periods": local_periods, **data})
            if len(tasks) == args.tasks:
                break
    assert len(tasks) == args.tasks
    k = np.arange(DEGREE+1)
    phase = 2*np.pi*np.outer(np.arange(FREQUENCIES), k)/FREQUENCIES
    packet = {
        "degree": DEGREE, "L": FREQUENCIES, "seed": args.seed, "batch": 131072, "tasks": tasks,
        "DC": (np.cos(phase)/FREQUENCIES/(RADIUS**k)).tolist(),
        "DS": (np.sin(phase)/FREQUENCIES/(RADIUS**k)).tolist(),
        "h_reduction_periods": periods, "attempted_draws": draw+1, "rejected_draws": rejections,
        "scope": "Selected numerical targets for a gap at C7 with A1,B1 positive. "
                 "Prime-period failures exclude only that low-coefficient branch. "
                 "Root preparation and mask ranking supply numerical exploration.",
    }
    with args.out.open("x") as output:
        json.dump(packet, output, sort_keys=True, separators=(",",":"))
        output.write("\n")
    print(json.dumps({"prepared": len(tasks), "draws": draw+1, "masks": sum(t["masks"] for t in tasks),
                      "rejections": dict(Counter(x["kind"] for x in rejections))}))


if __name__ == "__main__":
    main()
