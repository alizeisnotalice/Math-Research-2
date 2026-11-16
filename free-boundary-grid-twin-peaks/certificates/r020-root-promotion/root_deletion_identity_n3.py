#!/usr/bin/env python3
"""Full 3x3 permutation audit for the one-root deletion fibre identity."""
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


def common_neighbours(peaks):
    ans = set(range(N))
    for v in peaks:
        ans &= adj[v]
    return ans


def pair_type(peaks):
    a, b = peaks
    dx, dy = abs(a % n - b % n), abs(a // n - b // n)
    assert dx + dy == 2
    return 'd2_straight' if dx == 0 or dy == 0 else 'd2_diagonal'


A1 = 0
lhs = 0
rhs1 = rhs2_straight = rhs2_diagonal = rhs3 = rhs4 = rhs_other = 0
by_k = {}
for word in permutations(range(N)):  # highest label to lowest
    earlier, peaks = set(), []
    for v in word:
        if not (adj[v] & earlier):
            peaks.append(v)
        earlier.add(v)
    k = len(peaks)
    by_k[k] = by_k.get(k, 0) + 1
    if k == 1:
        A1 += 1
        rhs1 += 1 + len(adj[peaks[0]])
    elif k == 2:
        contribution = len(common_neighbours(peaks))
        if contribution:
            if pair_type(peaks) == 'd2_straight':
                rhs2_straight += contribution
            else:
                rhs2_diagonal += contribution
    elif k == 3:
        rhs3 += len(common_neighbours(peaks))
    elif k == 4:
        rhs4 += len(common_neighbours(peaks))
    else:
        rhs_other += len(common_neighbours(peaks))

lhs = N * A1
rhs = rhs1 + rhs2_straight + rhs2_diagonal + rhs3 + rhs4 + rhs_other
assert lhs == rhs
assert rhs_other == 0
print(json.dumps({
    'status': 'passed', 'n': n, 'A1': A1, 'N_A1': lhs,
    'peak_count_histogram': by_k,
    'terms': {'K1_root_and_neighbours': rhs1,
              'K2_distance2_straight': rhs2_straight,
              'K2_distance2_diagonal': rhs2_diagonal,
              'K3_common_star': rhs3,
              'K4_common_star': rhs4,
              'K_ge5_common_neighbour': rhs_other},
    'RHS': rhs,
}, indent=2, sort_keys=True))
