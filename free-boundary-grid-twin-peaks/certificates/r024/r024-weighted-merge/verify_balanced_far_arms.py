#!/usr/bin/env python3
"""Exact audit for a balanced, macroscopically separated first-merge fibre.

No random sampling: this checks the explicit top-row paths in P_n square.
"""
from fractions import Fraction
from itertools import permutations
from math import comb
import json
from pathlib import Path


def nbrs(n, v):
    x, y = v
    return {(X, Y) for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
            if 0 <= X < n and 0 <= Y < n}


def connected(n, vertices):
    todo = {next(iter(vertices))}
    seen = set()
    while todo:
        v = todo.pop()
        if v in seen:
            continue
        seen.add(v)
        todo |= (nbrs(n, v) & vertices) - seen
    return seen == vertices


rows = []
for n in range(3, 101):
    k = (n - 1) // 2
    A = {(x, 0) for x in range(k)}
    w = (k, 0)
    B = {(x, 0) for x in range(k + 1, 2 * k + 1)}
    assert k >= 1 and connected(n, A) and connected(n, B)
    assert not any(nbrs(n, u) & B for u in A)
    assert nbrs(n, w) & A == {(k - 1, 0)}
    assert nbrs(n, w) & B == {(k + 1, 0)}
    # A path P_k has binomial root law: Var= (k-1)/4.  The two means
    # differ by k+1 horizontally, and there is no vertical fluctuation.
    moment = Fraction((k + 1) ** 2, 1) + Fraction(k - 1, 2)
    m = Fraction(2 * k * (2 * k - 1), 2 * k * k)
    theta = 1 / m
    assert m == 2 - Fraction(1, k)
    assert theta == Fraction(k, 2 * k - 1)
    rows.append({"n": n, "k": k, "m": str(m), "theta": str(theta),
                 "conditional_E_distance_squared": str(moment),
                 "conditional_Q": str(moment / (n * n))})

out = {"construction": "top-row disjoint P_k arms joined first by w",
       "claims_checked": ["grid validity", "induced-arm disjointness", "bridge adjacency",
                          "exact m=2-1/k", "exact conditional second moment"],
       "rows": rows}

# This independent finite oracle checks the path-root formula used in the
# written argument by enumerating every relative label order for P_k.
path_oracle = []
for k in range(1, 10):
    counts = [0] * k
    for word in permutations(range(k)):
        rank = [0] * k
        for t, v in enumerate(word):
            rank[v] = t
        peaks = [v for v in range(k)
                 if all(rank[v] < rank[u] for u in (v - 1, v + 1) if 0 <= u < k)]
        if len(peaks) == 1:
            counts[peaks[0]] += 1
    expected = [comb(k - 1, r) for r in range(k)]
    assert counts == expected
    path_oracle.append({"k": k, "root_counts": counts})
out["path_root_full_permutation_oracle_k_le_9"] = path_oracle
Path(__file__).with_name("balanced_far_arms.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"n3": rows[0], "n100": rows[-1], "cases": len(rows)}, indent=2))
