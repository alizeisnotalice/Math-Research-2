#!/usr/bin/env python3
"""Exact finite checks for the degree-seven support and odd root-bound tables.

This script verifies only combinatorial/support and rational-polynomial
arithmetic used in round047's analytic proof.  It is not a counterexample
search and it makes no assertion about arbitrary cofactor degrees.
"""

from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json


DEGREE = 7
FORBIDDEN = {(i, i) for i in range(1, DEGREE)}
FORBIDDEN |= {(1, k) for k in (2, 3, 4)} | {(k, 1) for k in (2, 3, 4)}
FORBIDDEN |= {(2, 4), (4, 2)}


def conv(a, b):
    out = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def polynomial_boundary():
    p = [Fraction(1)]
    for _ in range(3):
        p = conv(p, [Fraction(1), Fraction(2)])
    for _ in range(3):
        p = conv(p, [Fraction(1), Fraction(1, 2)])
    p = conv(p, [Fraction(1), Fraction(1)])
    return p


def canonical(s):
    rev = tuple(sorted(DEGREE - i for i in s))
    return min(tuple(s), rev)


def main():
    all_sets = []
    kept = []
    for size in range(1, DEGREE):
        for s in combinations(range(1, DEGREE), size):
            r, t = min(s), DEGREE - max(s)
            row = {"S": list(s), "r": r, "t": t, "forbidden": (r, t) in FORBIDDEN}
            all_sets.append(row)
            if not row["forbidden"]:
                kept.append(row)
    orbits = {}
    for row in kept:
        key = canonical(tuple(row["S"]))
        orbits.setdefault(key, []).append(row["S"])
    boundary = polynomial_boundary()
    expected = [Fraction(1), Fraction(17, 2), Fraction(117, 4), Fraction(419, 8),
                Fraction(419, 8), Fraction(117, 4), Fraction(17, 2), Fraction(1)]
    assert boundary == expected
    expected_orbits = {(1,), (1, 2), (2,), (2, 4), (2, 3, 4), (3,)}
    assert set(orbits) == expected_orbits
    assert len(kept) == 12 and len(all_sets) == 63
    out = {
        "degree": DEGREE,
        "support_sets_total": len(all_sets),
        "forbidden_pairs": sorted(FORBIDDEN),
        "surviving_supports": kept,
        "reversal_orbits": {str(k): v for k, v in sorted(orbits.items())},
        "odd_boundary_polynomial": "(1+2x)^3(1+x/2)^3(1+x)",
        "odd_boundary_coefficients": [str(x) for x in boundary],
        "power_sum_upper_bound": "3*(2^s+2^(-s))+1 for every positive integer s",
        "scope": "Exact finite audit of the support classification and boundary polynomial only.",
    }
    target = Path(__file__).with_name("degree7_support_and_bounds_output.json")
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
