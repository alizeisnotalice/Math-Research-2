#!/usr/bin/env python3
"""Exact n=3..5 top-row first-merge weights using h(S) recursion."""
from functools import lru_cache
from fractions import Fraction
import json
from pathlib import Path

A2 = {
    3: 143808,
    4: 779299174080,
    5: 7063261481362958280000,
}


def one(n):
    N = n * n
    full = (1 << N) - 1
    adj = [0] * N
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                if 0 <= X < n and 0 <= Y < n:
                    adj[v] |= 1 << (Y * n + X)

    def boundary(s):
        z = 0
        q = s
        while q:
            bit = q & -q
            z |= adj[bit.bit_length() - 1]
            q -= bit
        return z & ~s

    @lru_cache(maxsize=None)
    def h(s):
        if s == full:
            return 1
        b = boundary(s)
        ans = 0
        while b:
            bit = b & -b
            ans += h(s | bit)
            b -= bit
        return ans

    top = (1 << n) - 1
    H = h(top)
    weights = []
    for j in range(1, n - 1):
        w = H * (1 << (n - 3)) * __import__("math").comb(n - 1, j)
        weights.append(w)
    F = sum(weights)
    formula = H * (1 << (n - 3)) * ((1 << (n - 1)) - 2)
    assert F == formula
    return {"n": n, "H_top_row": str(H), "weights_by_bridge_x": [str(x) for x in weights],
            "F": str(F), "A2": str(A2[n]), "F_over_A2": str(Fraction(F, A2[n])),
            "memoized_supersets": h.cache_info().currsize}


rows = [one(n) for n in (3, 4, 5)]
Path(__file__).with_name("top_row_arm_mass_n3_n5.json").write_text(json.dumps(rows, indent=2) + "\n")
print(json.dumps(rows, indent=2))
