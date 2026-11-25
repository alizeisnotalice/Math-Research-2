"""Recompute the deterministic finite certificate, or compare to its saved JSON."""

import argparse
import json
from fractions import Fraction as F
from math import comb, factorial, lcm
from pathlib import Path

from approximants import rows
from intervals import certify, harmonic
from sondow import case


def run():
    signs = [
        certify(2, -1, {}),
        certify(1, -1, {}),
        certify(10**10, -5772156649, {}, limit=32),
    ]
    assert [v["decision"] for v in signs] == ["positive", "negative", None]
    # At 4 bits, rounding confuses an integer boundary.
    floor_test = certify(0, -69, {2: 100}, "floor", 4, 128)
    assert floor_test["bits"] > 4 and floor_test["decision"] == 0
    values = list(rows(12))
    for row in values:
        if row["family"] == "Wolfs":
            n = row["n"]
            scale = factorial(n) * lcm(*range(1, n + 1))
            weights = [
                F(comb(n, k) ** 2 * comb(n + k, k), factorial(k)) for k in range(n + 1)
            ]
            # Differentiate x^(n+k) log x by first differentiating in the exponent.
            numerator = sum(
                w * (4 * harmonic(k) - 2 * harmonic(n - k) - harmonic(n + k))
                for k, w in enumerate(weights)
            )
            assert (scale * numerator, scale * sum(weights)) == (row["P"], row["Q"])
    cases = [case(n) for n in range(1, 13)]
    assert all(
        r["scaled_threshold_difference"]["decision"] == "positive"
        and r["frac_minus_dI"]["decision"] == "positive"
        and r["excluded_divisors_of"] == r["D"]
        for r in cases
    )
    return {
        "sign_tests": signs,
        "floor_test": floor_test,
        "approximants": values,
        "sondow": cases,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, sort_keys=True)
    if args.check:
        if json.loads(data) != json.loads(args.check.read_text()):
            raise SystemExit("certificate mismatch")
        print("verified 26 approximants, 12 denominator tests, 4 interval cases")
    else:
        print(data)
