#!/usr/bin/env python3
"""Direct subset-DP check of the inert-root factorization for n=4,5."""
import json
from functools import cache


def counts(n):
    V = {(x, y) for x in range(n) for y in range(n)}
    S = {(2, 0), (0, 1), (1, 1), (2, 1)}
    S |= {(x, 0) for x in range(3, n)}
    S |= {(0, y) for y in range(2, n)}
    A = {(0, 0), (1, 0)}
    B = sorted(V - S - A)
    ix = {v: i for i, v in enumerate(B)}

    def N(v):
        x, y = v
        return {(x + dx, y + dy) for dx, dy in
                ((1, 0), (-1, 0), (0, 1), (0, -1)) if (x + dx, y + dy) in V}

    F = {v for v in B if N(v) & S}
    r = (1, 2)
    assert N(r) & set(B) <= F

    def run(vertices, roots):
        vertices = list(vertices)
        jx = {v: i for i, v in enumerate(vertices)}
        root_mask = sum(1 << jx[v] for v in roots)
        adj = [sum(1 << jx[w] for w in N(v) if w in jx) for v in vertices]
        full = (1 << len(vertices)) - 1

        @cache
        def h(mask):
            if mask == full:
                return 1
            avail = root_mask
            for i in range(len(vertices)):
                if mask >> i & 1:
                    avail |= adj[i]
            avail &= ~mask
            return sum(h(mask | (1 << i)) for i in range(len(vertices)) if avail >> i & 1)

        H = h(0)
        D = 0
        for i, x in enumerate(vertices):
            if not (root_mask >> i & 1):
                continue
            for j, y in enumerate(vertices):
                if root_mask >> j & 1:
                    continue
                if adj[i] >> j & 1:
                    D += h((1 << i) | (1 << j))
        return H, D

    H, D = run(B, F)
    B0 = [v for v in B if v != r]
    H0, D0 = run(B0, F - {r})
    b = len(B)
    assert H == b * H0 and D == (b - 2) * D0
    return {"n": n, "b": b, "H": H, "D": D, "H_core": H0, "D_core": D0}


print(json.dumps([counts(n) for n in (4, 5)], indent=2))
