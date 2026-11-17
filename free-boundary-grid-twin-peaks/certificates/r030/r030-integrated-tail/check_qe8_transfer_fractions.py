#!/usr/bin/env python3
"""Exact rational checks for the two integer-rounding inequalities in qe8-transfer.md."""
from fractions import Fraction
from math import isqrt
import json
from pathlib import Path


def floor_a_over_sqrt_x(a: int, x: Fraction) -> int:
    """floor(a/sqrt(x)) without floating point, for x>0."""
    lo, hi = 0, a * max(x.denominator, 1) + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid * mid * x <= a * a:
            lo = mid
        else:
            hi = mid
    return lo


def check_event(p: Fraction) -> int:
    k = floor_a_over_sqrt_x(4, p)
    assert k >= 4
    # Kp <= 4 sqrt(p), and 16/(K+2) < 4 sqrt(p), squared exactly.
    assert k * k * p <= 16
    assert p * (k + 2) * (k + 2) > 16
    return k


def check_moment(q: Fraction) -> int:
    k = floor_a_over_sqrt_x(8, 2 * q)
    assert k >= 4
    # Kq <= 4 sqrt(2q), and 32/(K+2) < 4 sqrt(2q), squared exactly.
    assert k * k * q <= 32
    assert q * (k + 2) * (k + 2) > 32
    return k


event_inputs = [Fraction(a, b) for b in range(1, 65) for a in range(1, b + 1)]
moment_inputs = [Fraction(a, b) for b in range(1, 65) for a in range(1, 2 * b + 1)]
event_ks = [check_event(p) for p in event_inputs]
moment_ks = [check_moment(q) for q in moment_inputs]

out = {
    "method": "exact Fraction arithmetic; no floating point",
    "event_domain": "p=a/b, 1<=a<=b<=64",
    "moment_domain": "q=a/b, 1<=a<=2b, b<=64",
    "event_cases": len(event_inputs),
    "moment_cases": len(moment_inputs),
    "event_K_range": [min(event_ks), max(event_ks)],
    "moment_K_range": [min(moment_ks), max(moment_ks)],
    "checks": [
        "K=floor(4/sqrt(p)): K^2 p<=16<p(K+2)^2",
        "K=floor(8/sqrt(2q)): K^2 q<=32<q(K+2)^2",
    ],
}
target = Path(__file__).with_name("qe8-transfer-fraction-check.json")
target.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(out, indent=2, sort_keys=True))
