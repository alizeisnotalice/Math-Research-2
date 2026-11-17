#!/usr/bin/env python3
"""Complete exact N=7 bipartite check of R45=h3-3h4+2h5.

Every bipartite seven-vertex graph is isomorphic to one having a fixed
bipartition of size 1+6, 2+5, or 3+4.  All edge subsets of these complete
bipartite hosts are scanned; computations use only Python integers.
"""
from itertools import combinations, permutations
import json
from pathlib import Path

N = 7
WORDS = list(permutations(range(N)))


def hist(edges):
    adj = [0] * N
    for x, y in edges:
        adj[x] |= 1 << y; adj[y] |= 1 << x
    h = [0] * N
    for word in WORDS:
        r, p = word[0], word[1]
        if adj[p] & (1 << r):
            continue
        seen = (1 << r) | (1 << p)
        peaks = 2
        ell = 0
        for pos in range(2, N):
            v = word[pos]
            if not (adj[v] & seen):
                peaks += 1
                if peaks > 2:
                    break
            if not ell and (adj[p] & (1 << v)):
                ell = pos - 1  # L=J_full-2, positions are zero-based.
            seen |= 1 << v
        if peaks == 2 and ell:
            h[ell] += 1
    return h


families = []
checked = 0
bad = []
for left_size in (1, 2, 3):
    pairs = [(x, y) for x in range(left_size) for y in range(left_size, N)]
    nonzero = 0
    maxabs = 0
    for mask in range(1 << len(pairs)):
        es = [pair for k, pair in enumerate(pairs) if mask >> k & 1]
        h = hist(es)
        rel = h[3] - 3 * h[4] + 2 * h[5]
        checked += 1
        if rel:
            nonzero += 1
            if len(bad) < 10:
                bad.append({"left_size": left_size, "edges": es,
                            "hist": {str(i): h[i] for i in range(1, N) if h[i]},
                            "R45": rel})
        maxabs = max(maxabs, abs(rel))
    families.append({"fixed_bipartition": [left_size, N-left_size],
                     "edge_subsets": 1 << len(pairs),
                     "nonzero_R45": nonzero, "max_abs_R45": maxabs})

assert not bad
out = {"method": "complete exact permutation enumeration", "N": N,
       "definition": "R45=h3-3h4+2h5", "families": families,
       "graphs_scanned_up_to_relabelling": checked, "counterexamples": bad}
Path(__file__).with_name("exhaust-bipartite7-r45.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, indent=2, sort_keys=True))
