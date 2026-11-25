"""Exact Aptekarev and mixed-orthogonality approximants at x=1."""

from fractions import Fraction as F
from math import comb, factorial, gcd, lcm

from intervals import certify
from intervals import harmonic as H


def exact(x):
    if F(x).denominator != 1:
        raise ArithmeticError("integer required")
    return int(x)


def aptekarev(n):
    p = q = 0
    for k in range(n + 1):
        w = comb(n, k) ** 2 * factorial(n + k)
        q += w
        p += w * (H(n + k) + 2 * H(n - k) - 2 * H(k))
    return exact(p), q


def recurrence(initial, last):
    seq = list(initial)
    for n in range(2, last):
        a = (128 * n**3 + 40 * n * n - 82 * n - 45) * seq[n]
        b = n * n * (256 * n**3 - 240 * n * n + 64 * n - 7) * seq[n - 1]
        c = n * n * (n - 1) ** 2 * (16 * n + 1) * seq[n - 2]
        seq.append(exact(F(a - b + c, 16 * n - 15)))
    return seq[: last + 1]


def wolfs(n):
    p = q = 0
    for k in range(n + 1):
        w = F(comb(n, k) ** 2, factorial(k))
        v = comb(n + k, k)
        q += w * v
        p += w * (
            v * (3 * H(k) - 2 * H(n - k))
            + sum(F((-1) ** j * comb(n + k, n - j), j) for j in range(1, n + 1))
        )
    scale = factorial(n) * lcm(*range(1, n + 1))
    return exact(scale * p), exact(scale * q)


def rows(last):
    ps, qs = recurrence([0, 2, 31], last), recurrence([1, 3, 50], last)
    for n in range(last + 1):
        pair = aptekarev(n)
        if pair != (ps[n], qs[n]):
            raise ArithmeticError(f"recurrence mismatch at {n}")
        for name, (P, Q) in [("Aptekarev", pair), ("Wolfs", wolfs(n))]:
            g = gcd(P, Q)
            p, q = P // g, Q // g
            bound = certify(
                q,
                -p,
                {},
                start=max(64, q.bit_length() + 48),
                limit=max(512, q.bit_length() + 48),
            )
            lo, hi = F(bound["lo"]), F(bound["hi"])
            yield {
                "family": name,
                "n": n,
                "P": P,
                "Q": Q,
                "gcd": g,
                "p": p,
                "q": q,
                "form": bound,
                "error_interval": [str(lo / q), str(hi / q)],
            }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("last", type=int)
    last = parser.parse_args().last
    if last < 0:
        parser.error("nonnegative final index required")
    print(json.dumps(list(rows(last)), indent=2))
