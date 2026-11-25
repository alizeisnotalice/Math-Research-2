"""Sondow's integer exponents, exact floor and finite denominator tests."""

from collections import Counter
from math import comb, lcm

from approximants import exact
from intervals import certify
from intervals import harmonic as H


def factors(x):
    p = 2
    while p * p <= x:
        while x % p == 0:
            yield p
            x //= p
        p += 1
    if x > 1:
        yield x


def case(n):
    if n < 1:
        raise ValueError("n must be positive")
    d = lcm(*range(1, 2 * n + 1))
    D = d * comb(2 * n, n)
    A = exact(d * sum(comb(n, i) ** 2 * H(n + i) for i in range(n + 1)))
    logs = Counter()
    for k in range(1, n + 1):
        exponent = exact(
            2
            * d
            * sum(
                comb(n, i) ** 2 * (H(n - i) - H(i))
                for i in range(min(k - 1, n - k) + 1)
            )
        )
        for p in factors(n + k):
            logs[p] += exponent
    precision = max(64, sum(logs.values()).bit_length() + 48)
    fl = certify(0, 0, logs, "floor", precision, max(512, precision))
    j = fl["decision"]
    if j is None:
        raise ArithmeticError("fractional part unresolved")
    threshold = certify(0, -(2**n * j + 1), {p: 2**n * v for p, v in logs.items()})
    # frac(log S) - dI = A - floor(log S) - D*gamma.
    diff = certify(-D, A - j, {})
    if diff["decision"] is None:
        excluded = None
    else:
        excluded = D
    return {
        "n": n,
        "d": d,
        "D": D,
        "dA": A,
        "prime_exponents": dict(logs),
        "log_floor": fl,
        "scaled_threshold_difference": threshold,
        "frac_minus_dI": diff,
        "dI": certify(D, -A, logs),
        "excluded_divisors_of": excluded,
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    args = parser.parse_args()
    if args.n < 1:
        parser.error("positive index required")
    print(json.dumps(case(args.n), indent=2))
