#!/usr/bin/env python3
"""Independent exhaustive verification of the beta=1 deletion/insertion identity on P3 square P3."""
import itertools
import json
from collections import Counter

n = 3
N = n * n
verts = tuple(range(N))
adj = {v: set() for v in verts}
for y in range(n):
    for x in range(n):
        v = y * n + x
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < n and 0 <= yy < n:
                adj[v].add(yy * n + xx)
dist = {(u, v): abs(u % n - v % n) + abs(u // n - v // n) for u in verts for v in verts}

def peaks(order, available):
    pos = {x: i for i, x in enumerate(order)}  # order is descending labels
    return [x for x in available if all(pos[x] < pos[y] for y in adj[x] if y in available)]

direct = Counter()
two_peak_total = 0
for order in itertools.permutations(verts):
    pk = peaks(order, verts)
    if len(pk) != 2:
        continue
    two_peak_total += 1
    u = order[0]
    v = pk[0] if pk[1] == u else pk[1]
    t = N - order.index(v)
    M = max(N - order.index(x) for x in adj[v])
    deleted = tuple(x for x in order if x != v)
    pk_deleted = peaks(deleted, set(verts) - {v})
    b = max((N - order.index(x) for x in pk_deleted if x != u), default=0)
    if M - b == t - 1:
        direct[dist[u, v]] += 1

rhs = Counter()
per_root = Counter()
for v in verts:
    remaining = tuple(x for x in verts if x != v)
    for order in itertools.permutations(remaining):
        pk = peaks(order, set(remaining))
        if len(pk) != 1:
            continue
        u = pk[0]
        if u in adj[v]:
            continue
        rhs[dist[u, v]] += 1
        per_root[(v, u)] += 1

result = {
    "n": n,
    "A2_direct": two_peak_total,
    "beta_one_distance_hist_direct": dict(sorted(direct.items())),
    "rhs_delete_one_peak_distance_hist": dict(sorted(rhs.items())),
    "equal_by_distance": direct == rhs,
    "beta_one_total": sum(direct.values()),
    "rhs_total": sum(rhs.values()),
    "nonzero_delete_root_pairs": len(per_root),
}
print(json.dumps(result, indent=2, sort_keys=True))
if direct != rhs:
    raise SystemExit("beta=1 identity failed")
