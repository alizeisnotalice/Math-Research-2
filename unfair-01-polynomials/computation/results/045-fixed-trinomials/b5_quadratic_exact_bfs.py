#!/usr/bin/env python3
"""Exact finite-degree probe for B=1+v*x^2+x^5, v=(3-sqrt(5))/2.

It is evidence for this one v and degree cap only.  It does not claim an
unbounded exclusion or treat arbitrary quadratic algebraic integers.
"""
from fractions import Fraction as F

ZERO = (F(0), F(0))
ONE = (F(1), F(0))
V = (F(3, 2), F(-1, 2))


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))
def mul(x, y): return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])


def sign(x):
    """Exact sign of p+q sqrt(5)."""
    p, q = x
    if q == 0:
        return (p > 0) - (p < 0)
    if p == 0:
        return (q > 0) - (q < 0)
    if (p > 0) == (q > 0):
        return 1 if p > 0 else -1
    comparison = p*p - 5*q*q
    # For opposite signs, p+q sqrt(5) has p's sign iff |p|>|q|sqrt(5).
    if comparison == 0:
        raise AssertionError("impossible nontrivial rational sqrt(5) equality")
    return (p > 0) - (p < 0) if comparison > 0 else (q > 0) - (q < 0)


def unit_interval(x):
    return sign(x) >= 0 and sign(sub(ONE, x)) >= 0


def bit(x):
    return x == ZERO or x == ONE


def trailing_ok(a):
    """a has leading coefficient 1; test C_{m+1},...,C_{m+5} with A higher terms zero."""
    m = len(a) - 1
    for k in range(m + 1, m + 6):
        value = ZERO
        if 0 <= k - 2 <= m:
            value = add(value, mul(V, a[k - 2]))
        if 0 <= k - 5 <= m:
            value = add(value, a[k - 5])
        if not bit(value):
            return False
    return True


def main():
    cap = 40
    states = {(ONE,): ()}
    candidates = []
    counts = []
    for k in range(1, cap + 1):
        nxt = {}
        for a, word in states.items():
            am2 = a[k - 2] if k >= 2 else ZERO
            am5 = a[k - 5] if k >= 5 else ZERO
            for c in (0, 1):
                ak = sub(sub((F(c), F(0)), mul(V, am2)), am5)
                if unit_interval(ak):
                    new = a + (ak,)
                    nxt.setdefault(new, word + (c,))
                    if ak == ONE and trailing_ok(new):
                        candidates.append((k, new, word + (c,)))
        states = nxt
        counts.append(len(states))
    assert not candidates, candidates
    print("PASS (finite evidence only): no exact A through degree", cap,
          "for v=(3-sqrt(5))/2")
    print("state counts by degree 1..40:", counts)


if __name__ == "__main__":
    main()
