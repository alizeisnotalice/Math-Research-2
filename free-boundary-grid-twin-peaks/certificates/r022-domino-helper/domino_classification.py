#!/usr/bin/env python3
"""Exact finite audit of non-star independent peak sets in closed dominoes.

For every horizontal/vertical domino in P_n square (n=3,...,7), enumerate
independent subsets of its closed-neighbourhood union of sizes 3 and 4,
discard those with an open common neighbour, and classify by number of
leaves on the two sides.  This is a geometry audit, not an enumeration of
labellings.
"""
from collections import Counter
from itertools import combinations
import json


def audit(n):
    vertices = [(x, y) for x in range(n) for y in range(n)]
    vset = set(vertices)
    nbr = {v: {u for u in ((v[0] - 1, v[1]), (v[0] + 1, v[1]),
                            (v[0], v[1] - 1), (v[0], v[1] + 1)) if u in vset}
           for v in vertices}
    closed = {v: nbr[v] | {v} for v in vertices}
    counts = Counter()
    bad = []
    records = 0
    for a in vertices:
        for b in nbr[a]:
            if a > b:  # unoriented domino; the classification is side-symmetric
                continue
            left = nbr[a] - {b}
            right = nbr[b] - {a}
            h = closed[a] | closed[b]
            for size in (3, 4):
                for p in combinations(sorted(h), size):
                    p = set(p)
                    if any(y in nbr[x] for x, y in combinations(p, 2)):
                        continue
                    common = set(vertices)
                    for x in p:
                        common &= nbr[x]
                    if common:
                        continue
                    records += 1
                    centre_count = int(a in p) + int(b in p)
                    side_counts = tuple(sorted((len(p & left), len(p & right))))
                    counts[(size, centre_count, side_counts)] += 1
                    # Assertions constituting the claimed local classification.
                    if centre_count or side_counts not in ({(1, 2)} if size == 3
                                                            else {(1, 3), (2, 2)}):
                        bad.append({'a': a, 'b': b, 'P': sorted(p),
                                    'size': size, 'centres': centre_count,
                                    'side_counts': side_counts})
    assert not bad
    return {'n': n, 'nonstar_independent_sets': records,
            'classes': {f'K{size};centres={centres};sides={sides}': value
                        for (size, centres, sides), value in sorted(counts.items())}}


print(json.dumps({'status': 'passed', 'audits': [audit(n) for n in range(3, 8)]},
                 indent=2, sort_keys=True))
