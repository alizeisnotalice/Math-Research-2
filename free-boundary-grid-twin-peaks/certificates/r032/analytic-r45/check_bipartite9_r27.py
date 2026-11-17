#!/usr/bin/env python3
"""Exact selected N=9 bipartite checks of R27=h2-5h4+5h6-2h7."""
from itertools import permutations
import json
from pathlib import Path


def histogram(n, edges):
    adj = [0] * n
    for x, y in edges:
        adj[x] |= 1 << y; adj[y] |= 1 << x
    h = [0] * n
    for word in permutations(range(n)):
        r, p = word[:2]
        if adj[p] & (1 << r):
            continue
        seen = (1 << r) | (1 << p); peaks = 2; ell = 0
        for pos in range(2, n):
            v = word[pos]
            if not adj[v] & seen:
                peaks += 1
                if peaks > 2: break
            if not ell and adj[p] & (1 << v):
                ell = pos - 1
            seen |= 1 << v
        if peaks == 2 and ell:
            h[ell] += 1
    return h


def path(n): return [(i, i + 1) for i in range(n - 1)]
base = [(x, y) for x in (0, 1) for y in (2, 3, 4)]
graphs = {
    "P9": (9, path(9)),
    "K23_pendant_path4": (9, base + [(2,5),(5,6),(6,7),(7,8)]),
    "K23_two_arms": (9, base + [(2,5),(5,6),(3,7),(7,8)]),
    "K23_star_extension": (9, base + [(2,5),(5,6),(5,7),(5,8)]),
    "P7_two_center_leaves": (9, path(7) + [(3,7),(3,8)]),
    "C6_path3": (9, path(6) + [(5,0),(0,6),(6,7),(7,8)]),
}
out = {}
for name, (n, edges) in graphs.items():
    h = histogram(n, edges)
    out[name] = {"edges": edges,
                 "hist": {str(i): h[i] for i in range(1, n) if h[i]},
                 "R27": h[2] - 5*h[4] + 5*h[6] - 2*h[7]}
print(json.dumps(out, indent=2, sort_keys=True))
Path(__file__).with_name("bipartite9-r27.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
