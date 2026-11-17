#!/usr/bin/env python3
"""Exact checks of B_q2=2h1-B and the triangle-free A1=4h1-2B identity."""
from itertools import permutations
import json
from pathlib import Path


def peaks(word, adj):
    pos = [0] * len(adj)
    for i, v in enumerate(word):
        pos[v] = i
    return [v for v in range(len(adj)) if all(pos[v] < pos[u] for u in adj[v])]


def audit(n, edges):
    adj = [set() for _ in range(n)]
    for x, y in edges:
        adj[x].add(y); adj[y].add(x)
    a1 = b = h1 = bq2 = 0
    for word in permutations(range(n)):
        pp = peaks(word, adj)
        if len(pp) == 1:
            a1 += 1
        if len(pp) != 2 or set(word[:2]) != set(pp):
            continue
        b += 1
        r, p, x = word[:3]
        if p in adj[x]:
            h1 += 1
        if x in adj[r] and x in adj[p]:
            bq2 += 1
    triangle_sum = 0
    # h(S) = no-new-birth suffix count after all three vertices in S appeared.
    # Here it is used only to check the all-graph correction identity.
    from itertools import combinations
    for s in combinations(range(n), 3):
        if not (s[1] in adj[s[0]] and s[2] in adj[s[0]] and s[2] in adj[s[1]]):
            continue
        rem = [v for v in range(n) if v not in s]
        hs = 0
        for tail in permutations(rem):
            seen = set(s)
            ok = True
            for v in tail:
                if not (adj[v] & seen):
                    ok = False; break
                seen.add(v)
            hs += ok
        triangle_sum += hs
    assert bq2 == 2 * h1 - b
    assert a1 - 6 * triangle_sum == 4 * h1 - 2 * b
    return {"A1": a1, "B": b, "h1": h1, "B_q2": bq2,
            "triangle_suffix_sum": triangle_sum,
            "Bq2_equals_2h1_minus_B": bq2 == 2 * h1 - b,
            "general_corrected_identity": a1 - 6 * triangle_sum == 4 * h1 - 2 * b}


graphs = {
    "P5": (5, [(0,1),(1,2),(2,3),(3,4)]),
    "C5": (5, [(0,1),(1,2),(2,3),(3,4),(4,0)]),
    "K3": (3, [(0,1),(1,2),(0,2)]),
    "triangle_with_leaf": (4, [(0,1),(1,2),(0,2),(0,3)]),
    "grid_2x3": (6, [(0,1),(1,2),(3,4),(4,5),(0,3),(1,4),(2,5)]),
}
out = {name: audit(n, edges) for name, (n, edges) in graphs.items()}
out["method"] = "complete permutation enumeration with exact integer checks"
Path(__file__).with_name("rank3-identity-small-graphs.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(out, indent=2, sort_keys=True))
