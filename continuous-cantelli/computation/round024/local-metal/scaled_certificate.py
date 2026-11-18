"""Rational lower certificate for the scaled symmetric quartic integral."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json


lam = Q(7, 20)**4
stationary = Q(7, 20)
panels = 4096
sqrt_scale = 2**26
sum_scale = 10**12
normalization = Q(7977, 10000)
assert normalization**2 * Q(22, 7) < 2


def reciprocal_root_lower(s):
    """Bound 1/sqrt(s) below by rounding sqrt(s) upwards."""
    t = s.numerator * sqrt_scale**2
    d = s.denominator
    q = isqrt(t // d)
    if q*q*d < t:
        q += 1
    assert q*q*d >= t
    return Q(sqrt_scale, q)


def exp_lower(exponent):
    term = total = Q(1)
    for j in range(1, 16):
        term *= -exponent/j
        total += term
    return total


def lower_panel(a, b, tail):
    if tail:
        variance = 1 + lam*b**4
        exponent = b*b / (2*variance)
    else:
        variance = lam + b**4
        t = min(max(stationary, a), b)
        exponent = t*t/(2*(lam+t**4))
    lower = exp_lower(exponent) * reciprocal_root_lower(variance)
    assert lower > 0
    scaled = lower*sum_scale
    return scaled.numerator//scaled.denominator


def main():
    sums = [
        sum(lower_panel(Q(j, panels), Q(j+1, panels), tail) for j in range(panels))
        for tail in (False, True)
    ]
    bound = normalization * Q(sum(sums), panels * sum_scale)
    assert bound > Q(21, 10)
    summary = {
        "lambda": str(lam),
        "panels_per_half": panels,
        "integer_panel_sums": sums,
        "normalization_lower_bound": str(normalization),
        "lower_bound_exact": str(bound),
        "lower_bound_decimal": float(bound),
        "claimed_strict_lower_bound": "21/10",
    }
    Path("scaled-certificate.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
