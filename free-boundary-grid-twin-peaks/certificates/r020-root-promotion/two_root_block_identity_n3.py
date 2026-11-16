#!/usr/bin/env python3
"""Full 3x3 audit of the ordered-edge two-root promotion identity."""
from itertools import permutations
import json

n = 3
N = n * n
adj = []
for v in range(N):
    x, y = v % n, v // n
    adj.append({yy * n + xx for xx, yy in
                ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
                if 0 <= xx < n and 0 <= yy < n})
closed = [a | {v} for v, a in enumerate(adj)]
ordered_edges = [(a, b) for a in range(N) for b in adj[a]]


def cover_count(peaks):
    P = set(peaks)
    return sum(P <= (closed[a] | closed[b]) for a, b in ordered_edges)


def pair_distance(peaks):
    a, b = peaks
    return abs(a % n - b % n) + abs(a // n - b // n)


A1 = 0
by_k = {}
weighted_by_k = {}
k2_by_distance = {}
rhs = 0
for word in permutations(range(N)):
    earlier, peaks = set(), []
    for v in word:
        if not (adj[v] & earlier):
            peaks.append(v)
        earlier.add(v)
    k = len(peaks)
    e = cover_count(peaks)
    by_k[k] = by_k.get(k, 0) + 1
    weighted_by_k[k] = weighted_by_k.get(k, 0) + e
    if k == 2:
        d = pair_distance(peaks)
        k2_by_distance[d] = k2_by_distance.get(d, 0) + e
    if k == 1:
        A1 += 1
    rhs += e

lhs = N * (N - 1) * A1
assert lhs == rhs
assert all(value == 0 for d, value in k2_by_distance.items() if d > 3)
assert all(value == 0 for k, value in weighted_by_k.items() if k > 4)
print(json.dumps({
    'status': 'passed', 'n': n, 'A1': A1, 'N_falling_2_A1': lhs,
    'ordered_edges': len(ordered_edges), 'RHS': rhs,
    'peak_count_histogram': by_k,
    'weighted_by_peak_count': weighted_by_k,
    'K2_weighted_by_manhattan_distance': k2_by_distance,
}, indent=2, sort_keys=True))
