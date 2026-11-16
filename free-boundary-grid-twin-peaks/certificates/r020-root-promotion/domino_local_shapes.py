#!/usr/bin/env python3
"""Enumerate local independent peak sets in closed domino neighbourhoods."""
from itertools import combinations
import json

n = 5
V = [(x, y) for y in range(n) for x in range(n)]
N = {v: {w for w in V if abs(v[0]-w[0]) + abs(v[1]-w[1]) == 1} for v in V}


def independent(P):
    return all(q not in N[p] for p, q in combinations(P, 2))


def common(P):
    ans = set(V)
    for p in P:
        ans &= N[p]
    return ans


def d3pairs(P):
    return [pair for pair in combinations(P, 2)
            if abs(pair[0][0]-pair[1][0]) + abs(pair[0][1]-pair[1][1]) == 3]


def front_promotion_helper_for_K3(P):
    """Find p,z,q for the rigorous stable-front K3-to-K2 lemma."""
    if len(P) != 3:
        return None
    for p, q in d3pairs(P):
        r = next(iter(set(P) - {p, q}))
        for z in V:
            if z in P:
                continue
            if z in N[p] and z in N[r] and z not in N[q]:
                return {'p': p, 'z': z, 'q': q, 'deleted_peak': r}
            if z in N[q] and z in N[r] and z not in N[p]:
                return {'p': q, 'z': z, 'q': p, 'deleted_peak': r}
    return None


# Inspect a representative interior horizontal domino; all rotations/reflections
# of its two exceptional K4 patterns are geometrically equivalent.
a, b = (2, 2), (3, 2)
H = (N[a] | {a}) | (N[b] | {b})
records = []
for r in (2, 3, 4):
    for P in combinations(sorted(H), r):
        if not independent(P) or common(P):
            continue
        assert d3pairs(P)  # local geometry assertion used in r020 proof
        records.append({'P': P, 'size': r, 'd3_pairs': d3pairs(P),
                        'front_promotion_helper': front_promotion_helper_for_K3(P)})
assert all(x['front_promotion_helper'] is not None for x in records if x['size'] == 3)
print(json.dumps({'status': 'passed', 'domino': [a, b],
                  'nonstar_independent_records': records}, indent=2, sort_keys=True))
