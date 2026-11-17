#!/usr/bin/env python3
"""Exact small audit for the r037 fixed-gap fibre formula.

It is not a proof of the asymptotic theorem.  It checks finite indexing,
normalization, and a deliberately loose uniform O_t(1/q) bound.
"""

from fractions import Fraction
import platform
import sys


def falling(x: int, k: int) -> int:
    out = 1
    for s in range(k):
        out *= x - s
    return out


def p(a: int, b: int, t: int) -> Fraction:
    """Root-direction-averaged exact law of L=t in one merge fibre."""
    q = a + b
    assert a >= 1 and b >= 1 and t >= 1
    if q == 2:
        return Fraction(int(t == 1), 1)
    if a == 1 or b == 1:
        return Fraction(1, 2) if t in (1, q - 1) else Fraction(0, 1)
    m, x, y = q - 2, a - 1, b - 1
    if t > m:
        return Fraction(0, 1)
    return Fraction(x * falling(y, t - 1) + y * falling(x, t - 1),
                    2 * falling(m, t))


def g(theta: Fraction, t: int) -> Fraction:
    return (theta * (1 - theta) ** (t - 1)
            + (1 - theta) * theta ** (t - 1)) / 2


def main() -> None:
    normalization_cases = 0
    singleton_cases = 0
    approximation_cases = 0
    for q in range(2, 80):
        for a in range(1, q):
            b = q - a
            assert sum(p(a, b, t) for t in range(1, q)) == 1
            normalization_cases += 1
            if q > 2 and (a == 1 or b == 1):
                assert p(a, b, 1) == Fraction(1, 2)
                assert p(a, b, q - 1) == Fraction(1, 2)
                singleton_cases += 1
    for t in range(2, 10):
        for q in range(4 * t, 100):
            for a in range(1, q):
                b = q - a
                assert abs(p(a, b, t) - g(Fraction(a, q), t)) <= Fraction(32 * t * t, q)
                approximation_cases += 1
    print("status=PASS")
    print(f"python={sys.version.split()[0]}; implementation={platform.python_implementation()}")
    print(f"normalization_cases={normalization_cases}")
    print(f"singleton_cases={singleton_cases}")
    print(f"uniform_bound_cases={approximation_cases}")


if __name__ == "__main__":
    main()
