#!/usr/bin/env python3
"""Brute-force checks of equations (L) and (B) in leaf-bridge-identity.md."""
from itertools import permutations

def peak_set(word, adj):
    pos = {v: i for i, v in enumerate(word)}
    return {v for v in range(len(word)) if all(pos[v] < pos[w] for w in adj[v])}

def formula(n, edges, u):
    adj = [[] for _ in range(n)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    a2 = b = 0
    for word in permutations(range(n)):
        pset = peak_set(word, adj); k = len(pset); p = word.index(u) + 1; iu = u in pset
        r = word[0] in pset and word[1] in pset
        if k == 2: a2 += n + 1 if iu else n + 1 - p
        if k == 1 and not iu: a2 += p
        if k == 1 and not iu: b += 2
        if k == 2 and iu and p > 1: b += 2
        if k == 2 and r and iu: b += n - 1
        if k == 2 and r and not iu: b += n + 1 - p
    return a2, b

def direct(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    a2 = b = 0
    for word in permutations(range(n)):
        pset = peak_set(word, adj)
        if len(pset) == 2:
            a2 += 1
            b += int(word[0] in pset and word[1] in pset)
    return a2, b

def aggregate(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    a1 = a2 = b = marked = 0
    for word in permutations(range(n)):
        pset = peak_set(word, adj)
        if len(pset) == 1: a1 += 1
        if len(pset) == 2:
            a2 += 1
            marked += max(word.index(v) + 1 for v in pset)
            b += int(word[0] in pset and word[1] in pset)
    total = tuple(sum(x) for x in zip(*(direct(n + 1, edges + [(u, n)]) for u in range(n))))
    t = n * (n + 1) // 2
    predicted = ((t - 1) * a1 + (t + 1) * a2 + marked,
                 (t - 1) * b + 2 * a2 + 2 * (n - 1) * a1)
    assert total == predicted, (n, total, predicted)

cases = [
    (3, [(0, 1), (1, 2)], 2),
    (4, [(0, 1), (0, 2), (1, 3)], 3),
    (4, [(0, 1), (0, 2), (0, 3)], 1),
]
for n, e, u in cases:
    assert formula(n, e, u) == direct(n + 1, e + [(u, n)])
print('leaf-bridge formulas (L),(B): 3 cases passed')
aggregate(2, [(0, 1)])
aggregate(3, [(0, 1), (1, 2)])
print('all-attachment formulas (A),(C): edge and P3 passed')
