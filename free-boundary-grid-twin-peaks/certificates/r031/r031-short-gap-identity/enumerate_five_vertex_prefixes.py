#!/usr/bin/env python3
"""Exhaust all labelled 5-vertex graphs for the local h2-2h3 coefficient.

For an ordered five-prefix (r,p,x,y,z), r,p are the two births iff r is
first and rp is absent.  The last three entries must each have an earlier
neighbour.  A contributes to h2 when px is absent and py present; B
contributes to h3 when px,py are absent and pz present.  The coefficient is
1_A-2*1_B.  It is precisely the coefficient multiplying the no-new-birth
completion count f_G({r,p,x,y,z}).
"""
from itertools import combinations, permutations
from collections import defaultdict
import json
from pathlib import Path

N = 5
PAIRS = list(combinations(range(N), 2))
PERMS = list(permutations(range(N)))


def adjacent(mask, x, y):
    if x > y:
        x, y = y, x
    return bool(mask & (1 << PAIRS.index((x, y))))


def coefficient(mask):
    h2 = h3 = 0
    for r, p, x, y, z in PERMS:
        if adjacent(mask, r, p):
            continue
        prefix = [r, p]
        if not any(adjacent(mask, x, u) for u in prefix):
            continue
        prefix.append(x)
        if not any(adjacent(mask, y, u) for u in prefix):
            continue
        prefix.append(y)
        if not any(adjacent(mask, z, u) for u in prefix):
            continue
        if not adjacent(mask, p, x) and adjacent(mask, p, y):
            h2 += 1
        if not adjacent(mask, p, x) and not adjacent(mask, p, y) and adjacent(mask, p, z):
            h3 += 1
    return h2, h3


def is_bipartite(mask):
    colour = [None] * N
    for seed in range(N):
        if colour[seed] is not None:
            continue
        colour[seed] = 0
        todo = [seed]
        while todo:
            x = todo.pop()
            for y in range(N):
                if x == y or not adjacent(mask, x, y):
                    continue
                if colour[y] is None:
                    colour[y] = 1 - colour[x]
                    todo.append(y)
                elif colour[y] == colour[x]:
                    return False
    return True


def canonical(mask):
    # Exact brute-force isomorphism canonicalization; 5! * 1024 is tiny.
    bits = []
    for perm in PERMS:
        val = 0
        for k, (x, y) in enumerate(PAIRS):
            if adjacent(mask, perm[x], perm[y]):
                val |= 1 << k
        bits.append(val)
    return min(bits)


classes = {}
all_bipartite_zero = True
trianglefree_nonbipartite = []
for mask in range(1 << len(PAIRS)):
    h2, h3 = coefficient(mask)
    c = h2 - 2 * h3
    bip = is_bipartite(mask)
    if bip and c:
        all_bipartite_zero = False
    # Triangle test from triples.
    trianglefree = all(not (adjacent(mask, a, b) and adjacent(mask, a, c) and adjacent(mask, b, c))
                       for a, b, c in combinations(range(N), 3))
    if trianglefree and not bip:
        trianglefree_nonbipartite.append({"mask": mask, "c": c, "canonical": canonical(mask)})
    key = canonical(mask)
    if key not in classes:
        classes[key] = {"representative_mask": mask, "bipartite": bip,
                        "trianglefree": trianglefree, "h2": h2, "h3": h3,
                        "c": c, "labelled_members": 0}
    classes[key]["labelled_members"] += 1

def mask_of(edges):
    ans = 0
    for x, y in edges:
        if x > y:
            x, y = y, x
        ans |= 1 << PAIRS.index((x, y))
    return ans


c5 = mask_of([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)])
assert all_bipartite_zero
assert {x["canonical"] for x in trianglefree_nonbipartite} == {canonical(c5)}
assert {x["c"] for x in trianglefree_nonbipartite} == {20}

summary = {
    "method": "complete exact enumeration of all 2^10 labelled graphs and all 5! ordered prefixes",
    "definition": "c=h2_prefix-2*h3_prefix, with the no-new-birth completion factor excluded",
    "labelled_graphs": 1 << len(PAIRS),
    "isomorphism_classes": len(classes),
    "bipartite_all_c_zero": all_bipartite_zero,
    "trianglefree_nonbipartite_labelled_count": len(trianglefree_nonbipartite),
    "trianglefree_nonbipartite_canonical_classes": sorted({x["canonical"] for x in trianglefree_nonbipartite}),
    "trianglefree_nonbipartite_c_values": sorted({x["c"] for x in trianglefree_nonbipartite}),
    "classes": sorted(classes.values(), key=lambda x: x["representative_mask"]),
}
target = Path(__file__).with_name("five-vertex-prefix-enumeration.json")
target.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps({k: summary[k] for k in summary if k != "classes"}, indent=2, sort_keys=True))
