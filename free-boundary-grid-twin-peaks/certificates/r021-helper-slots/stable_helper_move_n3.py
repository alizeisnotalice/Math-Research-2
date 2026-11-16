#!/usr/bin/env python3
"""Full 3x3 audit of the stable helper-forward K3-to-K2 operation."""
from itertools import permutations
import json

n = 3
N = n * n
adj = []
for v in range(N):
    x, y = v % n, v // n
    adj.append({yy*n + xx for xx, yy in
                ((x-1,y), (x+1,y), (x,y-1), (x,y+1))
                if 0 <= xx < n and 0 <= yy < n})


def peaks(word):
    earlier, ans = set(), []
    for v in word:
        if not (adj[v] & earlier):
            ans.append(v)
        earlier.add(v)
    return ans


def choose_move(word, roots):
    pos = {v: i for i, v in enumerate(word)}
    rootset = set(roots)
    for p in roots:
        for q in roots:
            if p == q or abs(p % n-q % n)+abs(p//n-q//n) != 3:
                continue
            r = next(iter(rootset - {p, q}))
            if pos[p] >= pos[r]:
                continue
            for z in range(N):
                if z in rootset or pos[z] <= pos[r]:
                    continue
                if p in adj[z] and r in adj[z] and q not in adj[z]:
                    return p, z, q, r
    return None


total_k3 = moved = 0
fibre = {}
for word in permutations(range(N)):
    roots = peaks(word)
    if len(roots) != 3:
        continue
    total_k3 += 1
    move = choose_move(word, roots)
    if move is None:
        continue
    p, z, q, r = move
    old = word.index(z)
    before_r = word.index(r)
    new = list(word)
    new.pop(old)
    new.insert(before_r, z)  # z had been after r, so this is just before r.
    assert set(peaks(new)) == {p, q}
    moved += 1
    key = (tuple(new), p, z, q, r)
    fibre[key] = fibre.get(key, 0) + 1

print(json.dumps({
    'status': 'passed', 'n': n, 'K3_words': total_k3,
    'words_with_canonical_stable_helper_move': moved,
    'canonical_map_max_fibre_with_marked_p_z_q_r': max(fibre.values(), default=0),
    'canonical_map_fibre_histogram_with_marked_p_z_q_r': {
        str(k): sum(v == k for v in fibre.values()) for k in sorted(set(fibre.values()))},
    'canonical_map_distinct_marked_outputs': len(fibre),
}, indent=2))
