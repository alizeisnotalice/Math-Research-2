#!/usr/bin/env python3
"""Exact rational nullspaces among fixed-fibre orientation-averaged p_l."""
from fractions import Fraction
from math import prod, gcd
import json
from pathlib import Path


def falling(a, k):
    out = 1
    for j in range(k): out *= a-j
    return out


def p(x, y, ell):
    m = x + y
    return Fraction(y * falling(x, ell - 1) + x * falling(y, ell - 1), 2 * falling(m, ell))


def rref(a):
    a = [row[:] for row in a]
    nr, nc = len(a), len(a[0])
    pivots = []
    row = 0
    for col in range(nc):
        pivot = next((i for i in range(row, nr) if a[i][col]), None)
        if pivot is None: continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [v / scale for v in a[row]]
        for i in range(nr):
            if i != row and a[i][col]:
                f = a[i][col]
                a[i] = [u - f * v for u, v in zip(a[i], a[row])]
        pivots.append(col); row += 1
        if row == nr: break
    return a, pivots


def nullspace(a):
    rr, pivots = rref(a)
    nc = len(a[0]); frees = [j for j in range(nc) if j not in pivots]
    ans = []
    for free in frees:
        x = [Fraction(0) for _ in range(nc)]; x[free] = 1
        for i, pivot in reversed(list(enumerate(pivots))):
            x[pivot] = -sum(rr[i][j] * x[j] for j in frees)
        den = 1
        for z in x: den = den * z.denominator // gcd(den, z.denominator)
        ints = [int(z * den) for z in x]
        g = 0
        for z in ints: g = gcd(g, abs(z))
        ints = [z // g for z in ints]
        first = next(z for z in ints if z)
        if first < 0: ints = [-z for z in ints]
        ans.append(ints)
    return ans


out = {}
for upper in range(3, 13):
    # all sampled components exceed upper, so no small-support degeneracy enters.
    rows = [[p(x, y, ell) for ell in range(2, upper + 1)]
            for x in range(upper + 1, upper + 9)
            for y in range(upper + 1, upper + 9)]
    out[str(upper)] = {"ell_range": [2, upper], "nullspace": nullspace(rows)}
result = {"method": "exact Fraction row reduction on generic x,y samples", "results": out}
Path(__file__).with_name("generic-fibre-nullspace.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
