#!/usr/bin/env python3
"""Exhaustive 9! tau=2 ordered-word oracle for the L-distance observable."""
import itertools
import json

n = 3
N = n * n


def neighbours(v):
    x, y = v % n, v // n
    for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        xx, yy = x + dx, y + dy
        if 0 <= xx < n and 0 <= yy < n:
            yield yy * n + xx


adj = [set(neighbours(v)) for v in range(N)]


def md(a, b):
    return abs(a % n - b % n) + abs(a // n - b // n)


hist, sum_l, sum_l2, words = [0] * (2 * n - 1), 0, 0, 0
for word in itertools.permutations(range(N)):
    position = [0] * N
    for rank, vertex in enumerate(word):
        position[vertex] = rank
    births = [v for v in range(N) if not any(position[w] < position[v] for w in adj[v])]
    if len(births) != 2 or sorted(position[v] for v in births) != [0, 1]:
        continue
    vx, vy = word[0], word[1]
    assert vy not in adj[vx]
    best = next(rank + 1 for rank, vertex in enumerate(word) if rank >= 2 and vertex in adj[vy])
    L, d = best - 2, md(vx, vy)
    assert 1 <= L <= N - 3 and d >= 2
    words += 1
    hist[d] += L
    sum_l += L
    sum_l2 += L * L

assert words == 81984
assert hist == [0, 0, 63552, 64512, 15744]
assert sum_l == 143808
assert sum_l2 == 358080
print(json.dumps({"kind": "independent-n3-tau2-L-distance-oracle", "words": words,
                  "hist_sum_L": hist, "sum_L": sum_l, "sum_L2": sum_l2}, indent=2))
