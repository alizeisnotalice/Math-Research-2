#!/usr/bin/env python3
"""Construct legal stable-helper words with gain exactly 1/(N-2)."""
import json


def witness(n):
    assert n >= 5
    p, r, q, z = (1, 2), (2, 1), (4, 2), (2, 2)
    V = [(x, y) for y in range(n) for x in range(n)]
    neigh = {v: [(xx, yy) for xx, yy in
                 ((v[0]-1, v[1]), (v[0]+1, v[1]),
                  (v[0], v[1]-1), (v[0], v[1]+1))
             if 0 <= xx < n and 0 <= yy < n]
             for v in V}
    roots = [p, r, q]
    seen, queue, tail = set(roots), roots.copy(), []
    for v in queue:
        for u in neigh[v]:
            if u != z and u not in seen:
                seen.add(u); queue.append(u); tail.append(u)
    assert seen == set(V) - {z}
    word = roots + tail + [z]
    earlier, peaks = set(), []
    for v in word:
        if not (set(neigh[v]) & earlier):
            peaks.append(v)
        earlier.add(v)
    assert peaks == roots
    assert p in neigh[z] and r in neigh[z] and q not in neigh[z]
    P, R, Z = map(word.index, (p, r, z))
    assert (P, R, Z) == (0, 1, n*n-1)
    return {'n': n, 'roots': roots, 'helper': z,
            'positions': {'p': P, 'r': R, 'z': Z},
            'gain_numerator': R-P, 'gain_denominator': n*n-R-1}


print(json.dumps({'status': 'passed',
                  'witnesses': [witness(n) for n in range(5, 13)]}, indent=2))
