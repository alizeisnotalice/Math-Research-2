#!/usr/bin/env python3
"""Exact finite checks of the extended-L geometry; no random computation."""
from collections import deque
import json


def check(n):
    V = {(x, y) for x in range(n) for y in range(n)}
    S = {(2, 0), (0, 1), (1, 1), (2, 1)}
    S |= {(x, 0) for x in range(3, n)}
    S |= {(0, y) for y in range(2, n)}

    def N(A):
        return {(x + dx, y + dy) for x, y in A for dx, dy in
                ((1, 0), (-1, 0), (0, 1), (0, -1)) if (x + dx, y + dy) in V}

    q, seen = deque([next(iter(S))]), set()
    while q:
        x = q.popleft()
        if x in seen:
            continue
        seen.add(x)
        q.extend(N({x}) & S)
    assert seen == S
    degS = {x: len(N({x}) & S) for x in S}
    endpoints = {x for x, d in degS.items() if d == 1}
    assert set(degS.values()) <= {1, 2}
    assert sum(degS.values()) // 2 == len(S) - 1
    assert endpoints == {(n - 1, 0), (0, n - 1)}
    assert 2 * n - 2 == abs((n - 1) - 0) + abs(0 - (n - 1))

    R = V - S
    components = []
    while R:
        q, c = deque([next(iter(R))]), set()
        while q:
            x = q.popleft()
            if x in c:
                continue
            c.add(x)
            q.extend(N({x}) & R)
        R -= c
        components.append(c)
    components.sort(key=lambda c: (len(c), sorted(c)))
    A, B = components
    boundary = N(S) - S
    pairs = [(x, y) for x in boundary for y in N({x}) - S - N(S)]
    assert A == {(0, 0), (1, 0)}
    assert len(S) == 2 * n - 1
    assert len(B) == n * n - 2 * n - 1
    assert len(boundary) == 2 * n - 2
    assert len(boundary & B) == 2 * n - 4
    assert len(pairs) == 2 * n - 4 and all(x in B and y in B for x, y in pairs)
    return {"n": n, "size_S": len(S), "path_edges": len(S) - 1,
            "endpoints": sorted(endpoints), "m": n*n-len(S), "size_A": len(A),
            "size_B": len(B), "boundary": len(boundary),
            "B_boundary": len(boundary & B), "ordered_pairs": len(pairs)}


out = [check(n) for n in range(4, 21)]
print(json.dumps(out, indent=2))
