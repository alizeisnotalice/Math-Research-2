#!/usr/bin/env python3
"""Exact finite-fibre audit for r041; no asymptotic claim."""

from fractions import Fraction
import platform
import sys


def falling(x: int, k: int) -> int:
    value = 1
    for j in range(k):
        value *= x - j
    return value


def p(a: int, b: int, t: int) -> Fraction:
    q = a + b
    if q == 2:
        return Fraction(int(t == 1), 1)
    if a == 1 or b == 1:
        return Fraction(1, 2) if t in (1, q - 1) else Fraction()
    x, y, m = a - 1, b - 1, q - 2
    if t > m:
        return Fraction()
    return Fraction(x * falling(y, t - 1) + y * falling(x, t - 1),
                    2 * falling(m, t))


def main() -> None:
    cases = 0
    for q in range(2, 90):
        for a in range(1, q):
            b = q - a
            mean_from_law = sum(Fraction(t, 1) * p(a, b, t) for t in range(1, q))
            m = Fraction(q * (q - 1), 2 * a * b)
            theta, u = Fraction(a, q), Fraction(1, q)
            kernel = (1 - u) / (2 * theta * (1 - theta))
            assert sum(p(a, b, t) for t in range(1, q)) == 1
            assert mean_from_law == m == kernel
            cases += 1
    print("status=PASS")
    print(f"python={sys.version.split()[0]}; implementation={platform.python_implementation()}")
    print(f"fibres_checked={cases}")
    print("includes_q2=yes")


if __name__ == "__main__":
    main()
