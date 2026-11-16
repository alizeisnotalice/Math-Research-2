#!/usr/bin/env python3
"""Audit the exact local spin encoding for central-reflection label orbits."""
from itertools import combinations
import json
import random


def pairs(n):
    return [((x, y), (n - 1 - x, y))
            for y in range(n) for x in range(n // 2)]


def setup(n, label_pairs, spins):
    ps = pairs(n)
    pair_of, labels = {}, {}
    for i, (left, right) in enumerate(ps):
        lo, hi = sorted(label_pairs[i])
        high_left = bool((spins >> i) & 1)
        labels[left] = hi if high_left else lo
        labels[right] = lo if high_left else hi
        pair_of[left] = pair_of[right] = i
    neigh = {v: [(xx, yy) for xx, yy in
                 ((v[0] - 1, v[1]), (v[0] + 1, v[1]),
                  (v[0], v[1] - 1), (v[0], v[1] + 1))
                 if 0 <= xx < n and 0 <= yy < n]
             for v in labels}
    return ps, pair_of, labels, neigh


def direct_peaks(n, label_pairs, spins):
    _, _, labels, neigh = setup(n, label_pairs, spins)
    return {v for v in labels if all(labels[v] > labels[w] for w in neigh[v])}


def encoded_peaks(n, label_pairs, spins):
    ps, pair_of, labels, neigh = setup(n, label_pairs, spins)
    lows = [min(z) for z in label_pairs]
    highs = [max(z) for z in label_pairs]
    sig = [bool((spins >> i) & 1) for i in range(len(ps))]  # high on left
    ans = set()
    for u, (left, right) in enumerate(ps):
        for is_high in (True, False):
            side_left = sig[u] if is_high else not sig[u]
            vtx = left if side_left else right
            target = highs[u] if is_high else lows[u]
            ok = True
            for w in neigh[vtx]:
                v = pair_of[w]
                if v == u:
                    # Central-interface self pair: high beats low; low fails.
                    if not is_high:
                        ok = False
                    continue
                if lows[v] > target:
                    ok = False
                    break
                if lows[v] < target < highs[v]:
                    # The neighbour must take its low label on this side.
                    if is_high:
                        if sig[v] == sig[u]:
                            ok = False
                            break
                    else:
                        if sig[v] != sig[u]:
                            ok = False
                            break
            if ok:
                ans.add(vtx)
    return ans


def audit_orbit(n, label_pairs):
    for spins in range(1 << (n * n // 2)):
        assert direct_peaks(n, label_pairs, spins) == encoded_peaks(n, label_pairs, spins)


for a, b in combinations(range(1, 5), 2):
    audit_orbit(2, [(a, b), tuple(sorted(set(range(1, 5)) - {a, b}))])

witness = [(5, 12), (1, 7), (2, 11), (14, 15),
           (8, 10), (3, 16), (4, 9), (6, 13)]
audit_orbit(4, witness)
rng = random.Random(20261116)
for _ in range(100):
    vals = list(range(1, 17))
    rng.shuffle(vals)
    audit_orbit(4, [tuple(vals[2*i:2*i+2]) for i in range(8)])

print(json.dumps({'status': 'passed', 'n2_orbits': 6,
                  'n4_witness_orientations': 256,
                  'n4_sampled_orbits': 100}, indent=2))
