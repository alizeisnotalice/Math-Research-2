"""Certify the odd arc thresholds for the nine-body area and the almost-regular estimate."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

from flint import arb, ctx


def as_ball(q):
    v = Fraction(q)
    return arb(v.numerator) / v.denominator


def transferred(m):
    scale = 2 - 1 / (arb.pi() / (2 * m)).cos()
    return scale * scale * as_ball(Fraction(120697911165410052, 2**57))


def first_odd(target):
    # k parametrizes m=2k+1. The four finite inputs admit completions by 35 arcs.
    left, right = 17, 18
    threshold = as_ball(target)
    assert transferred(2 * left + 1) < threshold
    while not transferred(2 * right + 1) > threshold:
        right *= 2
    while right - left > 1:
        mid = (left + right) // 2
        value = transferred(2 * mid + 1)
        if value > threshold:
            right = mid
        else:
            assert value < threshold
            left = mid
    m = 2 * right + 1
    below, above = transferred(m - 2), transferred(m)
    assert below < threshold < above
    return dict(m=m, previous=below.str(44), certified=above.str(44))


def radius_stability(m):
    pi = arb.pi()
    optimal = 1 / (2 * (pi / (2 * m)).cos())
    radius = optimal + as_ball("0.00000001")
    gap = 2 * (1 / (2 * radius)).acos() - pi / m
    distance = ((m - 1) ** 2) * gap / 16
    assert gap > 0 and distance > as_ball("0.00000001")
    return {"arc_count": m, "radius_excess": "1/100000000", "upper_distance": distance.str(35)}


def check_dual(m):
    n = (m - 1) // 2
    checked = 0
    repeated = []
    for neutral in range(m):
        for positive in combinations([i for i in range(m) if i != neutral], n):
            t = [0 if i == neutral else 1 if i in positive else -1 for i in range(m)]
            values = [sum(t[i] * ((i - j - 1) % m) for i in range(m)) for j in range(m)]
            peak = max(values)
            assert peak <= n * (n + 1) and sum(values) == 0
            assert all((peak - x) % m == 0 for x in values)
            heights = [(peak - x) // m for x in values]
            assert sum(heights) == peak
            assert all(abs(heights[i] - heights[i-1]) <= 1 for i in range(m))
            if heights.count(0) > 1:
                assert peak <= n * n
                repeated.append(peak)
            checked += 1
    return dict(m=m, dual_vertices=checked, largest_peak_with_two_zeros=max(repeated))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    ctx.prec = 256
    targets = [".83750", ".837505", ".837509", ".8375099", ".83750999", ".83751"]
    results = {
        "arithmetic_bits": ctx.prec,
        "area_bound": str(Fraction(120697911165410052, 2**57)),
        "first_odd_for_this_transfer": {t: first_odd(t) for t in targets},
        "dual_checks": [check_dual(m) for m in (3, 5, 7, 9)],
        "stability_examples": [radius_stability(m) for m in (5, 7, 9)],
    }
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print([(t, item["m"]) for t, item in results["first_odd_for_this_transfer"].items()])
