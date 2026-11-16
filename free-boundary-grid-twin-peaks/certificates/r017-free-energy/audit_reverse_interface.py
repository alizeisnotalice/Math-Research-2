#!/usr/bin/env python3
"""Exact audit of the reverse right-interface root-forcing map on small grids."""

from collections import Counter
from itertools import permutations


def graph(w, h):
    V = tuple((x, y) for y in range(h) for x in range(w))
    adj = {v: {(v[0] + dx, v[1] + dy) for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))
               if (v[0] + dx, v[1] + dy) in V} for v in V}
    return V, adj


def peaks(word, adj):
    rank = {v: j for j, v in enumerate(word)}
    return {v for v in word if not any(rank[u] < rank[v] for u in adj[v])}


def D_order(w, h, q):
    I = tuple((w - 1, y) for y in range(h))  # fixed top-to-bottom connected order
    if q is None:
        return I
    x, y = q
    assert x < w - 1
    return I + tuple((xx, y) for xx in range(w - 2, x - 1, -1))


def audit(w, h):
    V, adj = graph(w, h)
    I = {(w - 1, y) for y in range(h)}
    fibres = Counter()
    source = 0
    for word in permutations(V):
        outside = peaks(word, adj) - I
        if len(outside) > 1:
            continue
        q = next(iter(outside), None)
        D = D_order(w, h, q)
        assert set(D) == (I if q is None else set(I) | {(x, q[1]) for x in range(q[0], w)})
        image = D + tuple(v for v in word if v not in set(D))
        assert peaks(image, adj) == {(w - 1, 0)}
        fibres[image] += 1
        source += 1
    m = w * h
    bound = (m + 1) * m ** (w + h - 1)
    assert max(fibres.values()) <= bound
    return {'rectangle': [w, h], 'source_T_right': source,
            'corner_onepeak_images': len(fibres), 'max_fibre': max(fibres.values()),
            'proved_fibre_bound': bound, 'all_images_corner_onepeak': True}


if __name__ == '__main__':
    print({'audits': [audit(3, 2), audit(3, 3)]})
