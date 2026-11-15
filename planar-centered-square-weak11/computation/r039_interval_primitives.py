#!/usr/bin/env python3
"""Exact dyadic root intervals and local span-rejection margin. Not an R39 certificate."""
from fractions import Fraction as F
from math import isqrt

def sqrt_interval(q, bits):
    q = F(q)
    if q < 0 or bits < 0:
        raise ValueError("q and bits must be nonnegative")
    scale = 1 << bits
    k = isqrt((q.numerator * scale * scale) // q.denominator)
    lo = F(k, scale)
    hi = lo if lo * lo == q else F(k + 1, scale)
    assert lo * lo <= q <= hi * hi
    return lo, hi

def rejection_radius(lam, mass, d0, u, H):
    lam, mass, d0, u, H = map(F, (lam, mass, d0, u, H))
    if lam <= 0 or mass < 0 or H <= 0:
        raise ValueError("invalid inputs")
    q = lam * d0 * d0 - mass
    if q <= 0:
        raise ValueError("strictly positive margin required")
    if u == 0:
        return H
    L = lam * abs(u) * (2 * abs(d0) + abs(u) * H)
    return min(H, q / (2 * L))

def selftest():
    cases = 0
    for n in range(0, 30):
        for d in range(1, 20):
            for bits in (0, 1, 4, 12):
                q = F(n, d)
                lo, hi = sqrt_interval(q, bits)
                assert lo * lo <= q <= hi * hi and hi - lo <= F(1, 1 << bits)
                cases += 1
    for d0 in (F(2), F(3), F(7, 3)):
        for u in (F(-1), F(0), F(1)):
            lam, mass, H = F(1), F(1), F(1, 2)
            delta = rejection_radius(lam, mass, d0, u, H)
            q = lam * d0 * d0 - mass
            for s in (-delta, F(0), delta):
                assert lam * (d0 + u * s)**2 - mass >= q / 2
    print(f"PASS root_enclosures={cases}, rejection_endpoint_checks=27; exact rational arithmetic")

if __name__ == "__main__":
    selftest()
