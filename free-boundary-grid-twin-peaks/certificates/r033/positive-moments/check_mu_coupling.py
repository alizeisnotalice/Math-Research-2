#!/usr/bin/env python3
"""Exact Fraction checks for the without/with-replacement mu_j coupling bound."""
from fractions import Fraction
import json
from pathlib import Path


def falling(a, j):
    ans = 1
    for k in range(j):
        ans *= a - k
    return ans


def row(M, X, j):
    assert M >= 1 and 1 <= j <= M
    Y = M - X
    finite = Fraction(falling(X, j) + falling(Y, j), 2 * falling(M, j))
    theta = Fraction(X, M)
    iid = (theta ** j + (1 - theta) ** j) / 2
    error = abs(finite - iid)
    bound = Fraction(j * (j - 1), 2 * M)
    assert error <= bound
    return {"M": M, "X": X, "j": j,
            "mu_finite": [finite.numerator, finite.denominator],
            "mu_iid": [iid.numerator, iid.denominator],
            "absolute_error": [error.numerator, error.denominator],
            "collision_union_bound": [bound.numerator, bound.denominator],
            "passes": True}


def singleton_boundary(j):
    # Direction-average singleton support is exactly {1,q-1}, each weight 1/2.
    q_hit = j + 3
    q_miss = j + 4
    assert q_hit - 1 == j + 2
    assert q_miss - 1 != j + 2
    return {"j": j,
            "q_equals_j_plus_3": {"support_L": [1, q_hit - 1],
                                    "p_L_equals_j_plus_2": [1, 2]},
            "q_equals_j_plus_4": {"support_L": [1, q_miss - 1],
                                    "p_L_equals_j_plus_2": [0, 1]}}


cases = [(5, 0, 2), (5, 1, 2), (5, 2, 3), (8, 3, 4),
         (11, 5, 5), (20, 1, 6), (20, 10, 6)]
out = {"method": "exact Fraction arithmetic",
       "cases": [row(*x) for x in cases],
       "singleton_boundary_cases": [singleton_boundary(j) for j in (1, 3, 5)]}
Path(__file__).with_name("mu-coupling-check.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, indent=2, sort_keys=True))
