#!/usr/bin/env python3
"""Exact brute force: B=exactly two peaks and its L=J-1 short-gap buckets.

Words are decreasing-label orders.  On B the first two entries are the two
peaks; L is the 0-based position, in the word after deleting the first entry,
of the first neighbour of the second entry.
"""
from itertools import permutations
from collections import Counter
import json
from pathlib import Path


def graph_path(n):
    return [(i, i + 1) for i in range(n - 1)]


def graph_cycle(n):
    return graph_path(n) + [(n - 1, 0)]


def graph_star(n):
    return [(0, i) for i in range(1, n)]


def graph_grid(rows, cols):
    ans = []
    for r in range(rows):
        for c in range(cols):
            v = r * cols + c
            if r + 1 < rows:
                ans.append((v, v + cols))
            if c + 1 < cols:
                ans.append((v, v + 1))
    return ans


def graph_t5():
    # Degree sequence (3,2,1,1,1).
    return [(0, 1), (0, 2), (0, 3), (1, 4)]


def graph_c4_leaf():
    return graph_cycle(4) + [(0, 4)]


def graph_k23():
    return [(x, y) for x in (0, 1) for y in (2, 3, 4)]


def graph_p3_union_k2():
    return [(0, 1), (1, 2), (3, 4)]


def count(n, edges):
    adj = [set() for _ in range(n)]
    for x, y in edges:
        adj[x].add(y)
        adj[y].add(x)
    hist = Counter()
    total = 0
    for word in permutations(range(n)):
        pos = [0] * n
        for i, x in enumerate(word):
            pos[x] = i
        peaks = [x for x in range(n) if all(pos[x] < pos[y] for y in adj[x])]
        if len(peaks) != 2 or set(word[:2]) != set(peaks):
            continue
        # In fact the set condition follows from K=2 plus top-two being peaks;
        # keeping it explicit makes the counted B slice unambiguous.
        p = word[1]
        first_neighbour_position = min(pos[x] for x in adj[p])
        L = first_neighbour_position - 1
        assert L >= 1
        hist[L] += 1
        total += 1
    return total, hist


families = {
    "P4": (4, graph_path(4)),
    "P5": (5, graph_path(5)),
    "P6": (6, graph_path(6)),
    "C4": (4, graph_cycle(4)),
    "C5": (5, graph_cycle(5)),
    "C6": (6, graph_cycle(6)),
    "S4": (4, graph_star(4)),
    "S5": (5, graph_star(5)),
    "S6": (6, graph_star(6)),
    "T5": (5, graph_t5()),
    "C4_leaf": (5, graph_c4_leaf()),
    "K2_3": (5, graph_k23()),
    "P3_union_K2": (5, graph_p3_union_k2()),
    "grid_2x2": (4, graph_grid(2, 2)),
    "grid_2x3": (6, graph_grid(2, 3)),
}
out = {"definition": "exact K=2, with the two first (highest) vertices the peaks"}
for name, (n, edges) in families.items():
    total, hist = count(n, edges)
    h2, h3 = hist[2], hist[3]
    out[name] = {"n": n, "edges": edges, "B": total,
                 "hist_L": {str(k): hist[k] for k in sorted(hist)},
                 "h2": h2, "h3": h3, "h2_minus_2h3": h2 - 2 * h3}
print(json.dumps(out, indent=2, sort_keys=True))
Path(__file__).with_name("small-graph-short-gap.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
