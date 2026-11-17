#!/usr/bin/env python3
"""Exact-Fraction audit for direction-averaged first-merge interleavings."""
from fractions import Fraction
from itertools import product
from math import prod
import sympy as sp


def falling(x: int, k: int) -> int:
    return prod(x - j for j in range(k))


def p(a: int, b: int, ell: int) -> Fraction:
    """P(L=ell), after averaging the choice of the lower root branch."""
    if a == 1 or b == 1:
        # For q>2: if the non-singleton branch is lower then L=1; if the
        # singleton is lower it receives its first neighbour only at merge,
        # so L=q-1.  The two root directions have equal weight.
        if a == b == 1:
            return Fraction(int(ell == 1))
        return Fraction(1, 2) if ell in (1, a + b - 1) else Fraction(0)
    x, y = a - 1, b - 1
    m = x + y
    if ell > m or ell < 1:
        return Fraction(0)
    return Fraction(falling(x, ell - 1) * y + falling(y, ell - 1) * x,
                    2 * falling(m, ell))


def main() -> None:
    # All branch pairs q <= 30; exact only.  q >= 7 should cancel.
    failures = []
    edge = []
    for a, b in product(range(1, 30), repeat=2):
        q = a + b
        c = p(a, b, 3) - 3 * p(a, b, 4) + 2 * p(a, b, 5)
        if q >= 7 and c != 0:
            failures.append((a, b, c))
        if q <= 6 and c:
            edge.append((a, b, q, p(a,b,3),p(a,b,4),p(a,b,5),c))
    assert not failures, failures[:5]
    assert edge == [
        (1, 3, 4, Fraction(1,2), Fraction(0), Fraction(0), Fraction(1,2)),
        (1, 4, 5, Fraction(0), Fraction(1,2), Fraction(0), Fraction(-3,2)),
        (1, 5, 6, Fraction(0), Fraction(0), Fraction(1,2), Fraction(1)),
        (2, 3, 5, Fraction(1,6), Fraction(0), Fraction(0), Fraction(1,6)),
        (2, 4, 6, Fraction(1,8), Fraction(1,8), Fraction(0), Fraction(-1,4)),
        (3, 1, 4, Fraction(1,2), Fraction(0), Fraction(0), Fraction(1,2)),
        (3, 2, 5, Fraction(1,6), Fraction(0), Fraction(0), Fraction(1,6)),
        (3, 3, 6, Fraction(1,6), Fraction(0), Fraction(0), Fraction(1,6)),
        (4, 1, 5, Fraction(0), Fraction(1,2), Fraction(0), Fraction(-3,2)),
        (4, 2, 6, Fraction(1,8), Fraction(1,8), Fraction(0), Fraction(-1,4)),
        (5, 1, 6, Fraction(0), Fraction(0), Fraction(1,2), Fraction(1)),
    ], edge
    for a in range(2, 35):
        assert p(1, a, 1) == Fraction(1,2)
        assert p(1, a, a) == Fraction(1,2)
        for t in range(2, 12):
            if t != a:
                assert p(1, a, t) == 0
    print("q>=7 cancellation: passed for 1<=a,b<30 (exact Fraction)")
    print("nonzero q<=6 cases:")
    for row in edge:
        print(row)

    # General odd-j relation from root-direction averaging.  This is an
    # exact finite rational check, not an asymptotic inference.
    for j in range(1, 10, 2):
        failures_j = []
        for a, b in product(range(1, 35), repeat=2):
            q = a + b
            if q >= j + 4:
                lhs = sum((-1)**i * Fraction(sp.binomial(j, i)) * p(a,b,i+2)
                          for i in range(j)) - 2*p(a,b,j+2)
                if lhs:
                    failures_j.append((a,b,lhs))
        assert not failures_j, (j, failures_j[:3])
    print("odd-j binomial relations: passed j=1,3,...,9, q>=j+4")

    # Symbolic verification of the new j=3 identity in x=a-1,y=b-1.
    x, y = sp.symbols("x y")
    m = x + y
    fall = lambda z, k: sp.prod(z-u for u in range(k))
    ps = lambda t: sp.factor((x*fall(y,t-1)+y*fall(x,t-1))/(2*fall(m,t)))
    assert sp.factor(ps(3)-3*ps(4)+2*ps(5)) == 0
    print("SymPy factor: p3-3*p4+2*p5 = 0")


if __name__ == "__main__":
    main()
