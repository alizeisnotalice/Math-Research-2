#!/usr/bin/env python3
"""Independent n=3 enumeration for the minimum-deletion identities."""
import itertools
import json
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("--n",type=int,choices=(2,3),default=3)
n=parser.parse_args().n
N = n * n
V = tuple(range(N))
adj = {v: set() for v in V}
for y in range(n):
    for x in range(n):
        v = y * n + x
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < n and 0 <= yy < n:
                adj[v].add(yy * n + xx)

def unique_peak(order, W):
    pos = {v: i for i, v in enumerate(order)}
    p = [v for v in W if all(pos[v] < pos[w] for w in adj[v] if w in W)]
    return p[0] if len(p) == 1 else None

a = {}
for v in V:
    W = tuple(x for x in V if x != v)
    for order in itertools.permutations(W):
        u = unique_peak(order, W)
        if u is not None:
            a[v, u] = a.get((v, u), 0) + 1

rhs_A1 = sum(a.values())
edge_root = sum(c for (v, u), c in a.items() if u in adj[v])
rhs_B1 = rhs_A1 - edge_root

A1 = 0
B1 = 0
for order in itertools.permutations(V):
    u = unique_peak(order, V)
    if u is None:
        continue
    A1 += 1
    v = order[-1]
    if u not in adj[v]:
        B1 += 1

out = {
    "n": n,
    "A1_direct": A1,
    "A1_punctured_sum": rhs_A1,
    "root_adjacent_minimum_count": edge_root,
    "B1_from_minimum_deletion": B1,
    "B1_punctured_nonadjacent_sum": rhs_B1,
    "passed": A1 == rhs_A1 and B1 == rhs_B1,
}
print(json.dumps(out, indent=2, sort_keys=True))
if not out["passed"]:
    raise SystemExit("identity failed")
