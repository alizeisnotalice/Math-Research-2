#!/usr/bin/env python3
"""Exact permutation checks of h3-3h4+2h5 on selected connected bipartite graphs."""
from itertools import permutations
from collections import Counter
import json
from pathlib import Path


def histogram(n, edges):
    adj = [0] * n
    for x, y in edges:
        adj[x] |= 1 << y; adj[y] |= 1 << x
    h = Counter(); b = 0
    for word in permutations(range(n)):
        seen = 0; peaks = 0; second = -1
        positions = [0] * n
        for i, v in enumerate(word):
            positions[v] = i
            if not (adj[v] & seen):
                peaks += 1
            seen |= 1 << v
        if peaks != 2 or not ((adj[word[0]] & 0) == 0):
            continue
        # first is automatically a peak; second must be a peak, exactly B.
        if adj[word[1]] & (1 << word[0]):
            continue
        b += 1
        p = word[1]
        j = min(positions[u] for u in range(n) if adj[p] >> u & 1)
        h[j - 1] += 1
    return b, h


def path(n): return [(i, i + 1) for i in range(n - 1)]
def add(n, base, es): return (n, base + es)

graphs = {
    "P5": (5, path(5)),
    "P6": (6, path(6)),
    "P7": (7, path(7)),
    "P6_leaf_at_2": add(7, path(6), [(2, 6)]),
    "P6_leaf_at_1": add(7, path(6), [(1, 6)]),
    "P5_two_leaves_center": add(7, path(5), [(2, 5), (2, 6)]),
    "C6_leaf": add(7, path(6) + [(5, 0)], [(0, 6)]),
    "K2_3_leaf": (6, [(x,y) for x in (0,1) for y in (2,3,4)] + [(2,5)]),
    "K2_3_two_leaves_same": (7, [(x,y) for x in (0,1) for y in (2,3,4)] + [(2,5),(2,6)]),
    "K2_3_two_leaves_split": (7, [(x,y) for x in (0,1) for y in (2,3,4)] + [(2,5),(3,6)]),
    "K2_4_leaf": (7, [(x,y) for x in (0,1) for y in (2,3,4,5)] + [(2,6)]),
    "K2_3_pendant_path3": (8, [(x,y) for x in (0,1) for y in (2,3,4)] + [(2,5),(5,6),(6,7)]),
    "K2_3_leaf_then_two": (8, [(x,y) for x in (0,1) for y in (2,3,4)] + [(2,5),(5,6),(5,7)]),
    "K2_3_path_at_core": (8, [(x,y) for x in (0,1) for y in (2,3,4)] + [(0,5),(5,6),(6,7)]),
}
out = {}
for name, (n, edges) in graphs.items():
    b, h = histogram(n, edges)
    rel45 = h[3] - 3*h[4] + 2*h[5]
    rel27 = h[2] - 5*h[4] + 5*h[6] - 2*h[7]
    out[name] = {"n":n, "edges":edges, "B":b,
                 "hist":{str(k): h[k] for k in sorted(h)},
                 "h3_minus_3h4_plus_2h5":rel45,
                 "h2_minus_5h4_plus_5h6_minus_2h7":rel27}
print(json.dumps(out, indent=2, sort_keys=True))
Path(__file__).with_name("bipartite-counterexample-search.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
