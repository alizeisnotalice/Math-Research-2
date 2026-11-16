#!/usr/bin/env python3
"""Cross-check small trees by brute force, then certify T_{3,3}."""
from __future__ import annotations

import hashlib
import itertools
import json
import pathlib
import time
from fractions import Fraction

import integer_tree_dp as dp


def brute(adj):
    n = len(adj)
    a1 = a2 = b = 0
    for word in itertools.permutations(range(n)):
        pos = {v: i for i, v in enumerate(word)}
        peaks = [v for v in range(n) if all(pos[v] < pos[w] for w in adj[v])]
        if len(peaks) == 1:
            a1 += 1
        if len(peaks) == 2:
            a2 += 1
            b += int(set(word[:2]) == set(peaks))
    return a1, a2, b


small = []
for n in range(3, 9):
    # path and a deterministic non-path tree where available
    path = [set() for _ in range(n)]
    for v in range(n - 1):
        path[v].add(v + 1); path[v + 1].add(v)
    cases = [path]
    if n >= 4:
        star = [set() for _ in range(n)]
        for v in range(1, n): star[0].add(v); star[v].add(0)
        cases.append(star)
    for adj in cases:
        got = (*dp.count_a1_a2(adj), dp.count_btop2(adj))
        want = brute(adj)
        assert got == want, (n, got, want)
        assert got[0] == dp.hook_a1(adj)
        small.append({'n': n, 'degrees': sorted(map(len, adj)), 'A1': got[0], 'A2': got[1], 'Btop2': got[2]})

g = dp.joined_binary(3, 3)
t0 = time.time()
a1, a2 = dp.count_a1_a2(g)
b = dp.count_btop2(g)
elapsed = time.time() - t0
expected_a2 = 4570935802581417984000000
expected_b = 1382099471396359372800000
assert (a2, b) == (expected_a2, expected_b), (a2, b)
ratio = Fraction(a2, b)
assert ratio > 3
out = {
    'status': 'passed',
    'method': 'independent exact integer rank-interleaving DP; fixed-top-two marked-forest integer no-peak DP',
    'source_sha256': hashlib.sha256(pathlib.Path('integer_tree_dp.py').read_bytes()).hexdigest(),
    'small_all_permutation_cases': small,
    'T_3_3': {'vertices': len(g), 'max_degree': max(map(len, g)), 'A1': str(a1), 'A2': str(a2), 'Btop2': str(b), 'ratio': str(ratio), 'ratio_float': float(ratio)},
    'wall_seconds': elapsed,
    'scope': 'Finite auxiliary counterexample to A2 <= Delta*Btop2 on all max-degree-3 trees; not a grid twin-peaks counterexample.'
}
pathlib.Path('receipt.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
