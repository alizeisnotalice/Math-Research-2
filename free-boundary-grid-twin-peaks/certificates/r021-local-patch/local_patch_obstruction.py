#!/usr/bin/env python3
"""Exact finite obstruction to an all-exterior-conditions local peak repair.

On P_7 square, fix the closed-neighbourhood patch of a horizontal domino.
External vertices form a rooted increasing forest into four prescribed
boundary threats.  Only the eight patch labels are permuted; all 8! cases
are exhaustively peak-tested.
"""
from collections import deque
from itertools import permutations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
n = int(sys.argv[1]) if len(sys.argv) == 2 else 7
if n < 7:
    raise SystemExit("this fixed local configuration requires n >= 7")
N = n * n


def vid(x, y):
    return y * n + x


def xy(v):
    return [v % n, v // n]


def neighbors(v):
    x, y = v % n, v // n
    return [vid(xx, yy) for xx, yy in ((x - 1, y), (x + 1, y),
            (x, y - 1), (x, y + 1)) if 0 <= xx < n and 0 <= yy < n]


# Domino a--b, its union of closed neighbourhoods, four pairwise
# nonadjacent leaves, and their unique patch neighbours x_i outside P.
a, b = vid(3, 3), vid(4, 3)
patch = sorted(set([a, b] + neighbors(a) + neighbors(b)))
leaves = [vid(2, 3), vid(3, 2), vid(5, 3), vid(4, 4)]
threats = [vid(1, 3), vid(3, 1), vid(6, 3), vid(4, 5)]
outside = sorted(set(range(N)) - set(patch))
assert len(patch) == 8 and len(outside) == N - 8
for leaf, threat in zip(leaves, threats):
    assert set(neighbors(threat)) & set(patch) == {leaf}
assert all(abs(xy(u)[0] - xy(v)[0]) + abs(xy(u)[1] - xy(v)[1]) > 1
           for i, u in enumerate(leaves) for v in leaves[i + 1:])


def fixed_exterior():
    """Assign exterior ranks: every non-threat has a higher exterior parent."""
    dist = {v: 0 for v in threats}
    parent = {}
    q = deque(threats)
    while q:
        u = q.popleft()
        for v in neighbors(u):
            if v in outside and v not in dist:
                dist[v] = dist[u] + 1
                parent[v] = u
                q.append(v)
    assert set(dist) == set(outside)
    lab = [None] * N
    # Larger values are closer to one of the four threats.
    nonroots = sorted((v for v in outside if v not in threats),
                      key=lambda v: (-dist[v], v))
    assert len(nonroots) == N - 12
    for rank, v in enumerate(nonroots, start=1):
        lab[v] = rank
    for rank, v in enumerate(threats, start=N - 7):
        lab[v] = rank
    assert all(lab[parent[v]] > lab[v] for v in nonroots)
    return lab, parent, dist


def peaks(labels):
    return tuple(v for v in range(N)
                 if all(labels[v] > labels[u] for u in neighbors(v)))


def common_open_neighbors(vertices):
    result = set(neighbors(vertices[0]))
    for v in vertices[1:]:
        result &= set(neighbors(v))
    return result


def ordered_dominating_edges(vertices):
    answer = []
    for u in range(N):
        closed_u = {u, *neighbors(u)}
        for v in neighbors(u):
            if all(p in closed_u | {v, *neighbors(v)} for p in vertices):
                answer.append((u, v))
    return answer


def main():
    labels, parent, dist = fixed_exterior()
    patch_labels = tuple(range(N - 11, N - 7)) + tuple(range(N - 3, N + 1))
    assert set(range(1, N + 1)) == set(x for x in labels if x is not None) | set(patch_labels)
    for threat in threats:
        assert all(labels[u] < labels[threat]
                   for u in neighbors(threat) if u in outside)
    source_value = {}
    for v, value in zip((v for v in patch if v not in leaves), range(N - 11, N - 7)):
        source_value[v] = value
    for v, value in zip(leaves, range(N - 3, N + 1)):
        source_value[v] = value
    source_order = tuple(source_value[v] for v in patch)
    source_labels = labels[:]
    for v, value in zip(patch, source_order):
        source_labels[v] = value
    source_peaks = peaks(source_labels)
    source_c = common_open_neighbors(source_peaks)
    source_edges = ordered_dominating_edges(source_peaks)
    assert set(source_peaks) == set(leaves)
    assert not source_c
    assert (a, b) in source_edges and (b, a) in source_edges
    min_k, examples, histogram = N + 1, [], {}
    # Stronger invariant: each non-threat exterior vertex is never a peak,
    # independently of the patch permutation.
    for order in permutations(patch_labels):
        trial = labels[:]
        for v, value in zip(patch, order):
            trial[v] = value
        p = peaks(trial)
        k = len(p)
        histogram[k] = histogram.get(k, 0) + 1
        if k < min_k:
            min_k, examples = k, [(order, p)]
        elif k == min_k:
            examples.append((order, p))
    assert min_k >= 3
    for order, p in examples:
        trial = labels[:]
        for v, value in zip(patch, order):
            trial[v] = value
        assert all(v not in p for v in outside if v not in threats)
    out = {
        "model": f"P_{n} square; exterior labels fixed; all 8! patch-label permutations tested",
        "domino": [xy(a), xy(b)], "patch": [xy(v) for v in patch],
        "selected_leaves": [xy(v) for v in leaves],
        "threats": [xy(v) for v in threats],
        "patch_label_multiset": list(patch_labels),
        "fixed_exterior_labels": {str(v): labels[v] for v in outside},
        "forest_parent": {str(v): parent[v] for v in parent},
        "forest_distance": {str(v): dist[v] for v in outside},
        "permutations": 40320, "minimum_peak_count": min_k,
        "peak_count_histogram": histogram,
        "local_domino_source": {
            "patch_values_in_sorted_patch_order": list(source_order),
            "peaks": [xy(v) for v in source_peaks],
            "peak_vertices": list(source_peaks),
            "peak_count": len(source_peaks),
            "independent": all(v not in neighbors(u) for i, u in enumerate(source_peaks)
                               for v in source_peaks[i + 1:]),
            "open_common_neighbors": [xy(v) for v in sorted(source_c)],
            "ordered_dominating_edges": [[xy(u), xy(v)] for u, v in source_edges],
        },
        "minimum_examples": [
            {"patch_values_in_sorted_patch_order": list(order),
             "peaks": [xy(v) for v in p], "peak_vertices": list(p)}
            for order, p in examples[:20]],
        "interpretation": "finite local obstruction only; not an original-model asymptotic counterexample",
    }
    (HERE / f"local-patch-obstruction-n{n}.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
