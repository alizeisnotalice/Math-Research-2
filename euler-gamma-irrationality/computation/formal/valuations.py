"""Prime-offset congruences; exact pairs and rational residual certificates."""

import argparse
import json
import sys
from fractions import Fraction
from math import comb, factorial, gcd, lcm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from approximants import aptekarev, recurrence, wolfs
from intervals import certify


def rivoal(n):
    initials = ([-1, 4, Fraction(77, 4)], [1, 7, Fraction(65, 2)])
    for v in initials:
        for j in range(n - 2):
            a, b, c = v[-3:]
            v.append(
                Fraction(
                    (j + 3) * (8 * j + 11) * (24 * j * j + 145 * j + 215) * c
                    - (8 * j + 27) * (24 * j**3 + 105 * j * j + 124 * j + 25) * b
                    + (j + 2) ** 2 * (8 * j + 19) * (8 * j + 27) * a,
                    (j + 3) ** 2 * (8 * j + 11) * (8 * j + 19),
                )
            )
    scale = factorial(n) * lcm(*range(1, n + 1))
    p, q = (Fraction(v[n]) * scale for v in initials)
    assert p.denominator == q.denominator == 1
    return int(p), int(q)


def run(route):
    ps = recurrence((0, 2, 31), 64)
    qs = recurrence((1, 3, 50), 64)
    primes = [p for p in range(7, 65) if all(p % i for i in range(2, int(p**0.5) + 1))]
    data = []
    for n in range(1, 65):
        P, Q = {"apt": aptekarev, "mixed": wolfs, "riv": rivoal}[route](n)
        if route == "apt":
            assert (P, Q) == (ps[n], qs[n])
        hits = []
        for p in primes:
            r = n - p
            if r < 0 or 2 * r + 1 >= p or route == "riv":
                continue
            if route == "apt":
                h = qs[r]
                assert Q % p == 0 and (Q // p + h) % p == 0 and (P + 3 * h) % p == 0
            else:
                h = sum(
                    comb(r, k) ** 2 * comb(r + k, k) * (factorial(r) // factorial(k))
                    for k in range(r + 1)
                )
                d = lcm(*range(1, n + 1))
                assert (
                    Q % d == 0
                    and (Q // d - 2 * h) % p == 0
                    and (P - 5 * (d // p) * h) % p == 0
                )
            hits.append([p, r, h % p])
        g = gcd(P, Q)
        q = Q // g
        bits = q.bit_length() + 100
        data.append(
            dict(
                n=n,
                P=P,
                Q=Q,
                gcd=g,
                p=P // g,
                q=q,
                tests=hits,
                form=certify(q, -(P // g), {}, start=bits, limit=max(bits, 512)),
            )
        )
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("route", choices=["apt", "mixed", "riv"])
    args = parser.parse_args()
    print(json.dumps(run(args.route), indent=2))
