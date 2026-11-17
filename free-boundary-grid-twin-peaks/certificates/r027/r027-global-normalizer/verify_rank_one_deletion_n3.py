#!/usr/bin/env python3
"""Exact rank-one deletion audit for P_n square, n=2,3."""
from itertools import permutations
from fractions import Fraction
import json
from pathlib import Path


def audit(n: int) -> dict:
    N = n * n
    adj = [set() for _ in range(N)]
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                if 0 <= X < n and 0 <= Y < n:
                    adj[v].add(Y * n + X)
    all_two = top_two = forward_c_total = forward_c_root_out = 0
    by_degree: dict[int, list[int]] = {}
    for word in permutations(range(N)):
        rank = [0] * N
        for i, v in enumerate(word):
            rank[v] = i
        peaks = {v for v in range(N) if all(rank[v] < rank[u] for u in adj[v])}
        if len(peaks) != 2:
            continue
        all_two += 1
        v = word[0]
        hpeaks = {
            z for z in range(N) if z != v
            and all(rank[z] < rank[u] for u in adj[z] if u != v)
        }
        assert len(hpeaks - adj[v]) == 1
        forward_c_total += 1
        root_out = word[1] not in adj[v]
        if root_out:
            top_two += 1
            forward_c_root_out += 1
        bucket = by_degree.setdefault(len(adj[v]), [0, 0])
        bucket[0] += 1
        bucket[1] += root_out
    # Independently enumerate each punctured graph rather than relying only on
    # deletion of already-known two-peak words.
    c_total = c_root_out = 0
    fibre_preimages: dict[tuple[int, tuple[int, ...]], set[int]] = {}
    for v in range(N):
        H = [z for z in range(N) if z != v]
        for word in permutations(H):
            rank = {z: i for i, z in enumerate(word)}
            hpeaks = {
                z for z in H
                if all(rank[z] < rank[u] for u in adj[z] if u != v)
            }
            if len(hpeaks - adj[v]) == 1:
                c_total += 1
                c_root_out += word[0] not in adj[v]
                p = next(iter(hpeaks - adj[v]))
                old_position = word.index(p) + 1
                promoted = (p,) + tuple(z for z in word if z != p)
                j = next(i + 1 for i, z in enumerate(promoted) if z in adj[p])
                assert old_position < j
                fibre_preimages.setdefault((v, promoted), set()).add(old_position)
    assert all_two == forward_c_total == c_total
    assert top_two == forward_c_root_out == c_root_out
    assert sum(len(x) for x in fibre_preimages.values()) == c_total
    for (v, word), positions in fibre_preimages.items():
        assert word[0] not in adj[v]
        j = next(i + 1 for i, z in enumerate(word) if z in adj[word[0]])
        assert positions == set(range(1, j))
    return {
        "n": n, "A2": all_two, "Btop2": top_two,
        "B_over_A2": str(Fraction(top_two, all_two)),
        "by_deleted_max_degree": {
            str(d): {"C": a, "root_out": b, "fraction": str(Fraction(b, a))}
            for d, (a, b) in sorted(by_degree.items())
        },
    }


out = {"kind": "rank-one-deletion-exhaustive-audit", "results": [audit(2), audit(3)]}
# A regression check for the failed root/outside-peak swap documented in the
# accompanying proof note.  Here v=0 and the two words are words of G-v.
adj3 = [set() for _ in range(9)]
for y in range(3):
    for x in range(3):
        z = y * 3 + x
        for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
            if 0 <= X < 3 and 0 <= Y < 3:
                adj3[z].add(Y * 3 + X)
def local_peaks(word):
    rank = {z: i for i, z in enumerate(word)}
    return {z for z in word if all(rank[z] < rank[u] for u in adj3[z] if u != 0)}
assert local_peaks((1, 2, 3, 4, 6, 8, 5, 7)) == {1, 3, 8}
assert local_peaks((8, 2, 3, 4, 6, 1, 5, 7)) == {2, 3, 8}
out["root_swap_counterexample_checked"] = True
Path(__file__).with_name("rank_one_deletion_n2_n3.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
