#!/usr/bin/env python3
"""Exact audit for the sharp vertical-interface family in r011.

For n=4,k=2 all binomial(16,8)=12870 relative-order shuffles are checked.
For n=5,6 every admissible cut is checked symbolically from the two side
orders: roots/safe-neighbourhoods, one internal higher neighbour per interface
vertex, and no other cross edges imply theta=2^{-n} for every shuffle.
"""

from fractions import Fraction
from itertools import combinations
from math import comb
import json
from pathlib import Path
import sys


def grid(n):
    return [tuple(
        yy * n + xx
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))
        if 0 <= (xx := x + dx) < n and 0 <= (yy := y + dy) < n
    ) for y in range(n) for x in range(n)]


def side_orders(n, k):
    # Lists are in strictly decreasing relative-label order.
    a = [x for x in range(k)]
    a += [y * n + (k - 1) for y in range(1, n)]
    a += [y * n + x for y in range(1, n) for x in range(k - 2, -1, -1)]

    b = [n - 1 - x for x in range(n - k)]
    b += [y * n + k for y in range(1, n)]
    b += [y * n + x for y in range(1, n) for x in range(k + 1, n)]
    assert set(a).isdisjoint(b) and set(a) | set(b) == set(range(n * n))
    return a, b


def labels_from_shuffle(a, b, a_positions):
    n_vertices = len(a) + len(b)
    # Scan global labels high-to-low: at an A-position consume next A order.
    ia = ib = 0
    label = [None] * n_vertices
    a_positions = set(a_positions)
    for high_rank in range(n_vertices):
        v = a[ia] if high_rank in a_positions else b[ib]
        ia += high_rank in a_positions
        ib += high_rank not in a_positions
        label[v] = n_vertices - 1 - high_rank
    assert ia == len(a) and ib == len(b)
    return label


def theta_and_roots(n, label):
    adj = grid(n)
    higher = [{u for u in adj[v] if label[u] > label[v]} for v in range(n * n)]
    roots = [v for v in range(n * n) if not higher[v]]
    return roots, higher


def structural_check(n, k):
    adj = grid(n)
    a, b = side_orders(n, k)
    aset, bset = set(a), set(b)
    arank = {v: i for i, v in enumerate(a)}
    brank = {v: i for i, v in enumerate(b)}
    p, q = a[0], b[0]
    assert all(u in aset for u in adj[p])
    assert all(u in bset for u in adj[q])
    for side, ranks, root in ((aset, arank, p), (bset, brank, q)):
        for v in side:
            internal_higher = [u for u in adj[v] if u in side and ranks[u] < ranks[v]]
            assert (v == root and not internal_higher) or (v != root and internal_higher)
    for y in range(n):
        av, bv = y * n + (k - 1), y * n + k
        assert sum(u in aset and arank[u] < arank[av] for u in adj[av]) == 1
        assert sum(u in bset and brank[u] < brank[bv] for u in adj[bv]) == 1
    return {"n": n, "k": k, "root_pair": [p, q], "root_distance": n - 1,
            "interface_edges": n, "forced_theta": f"1/{2 ** n}"}


def exhaustive_n4():
    n, k = 4, 2
    a, b = side_orders(n, k)
    checked = 0
    for apos in combinations(range(n * n), len(a)):
        label = labels_from_shuffle(a, b, apos)
        roots, higher = theta_and_roots(n, label)
        assert roots == [a[0], b[0]]
        theta = Fraction(1)
        for y in range(n):
            av, bv = y * n + k - 1, y * n + k
            lower = av if label[av] < label[bv] else bv
            assert len(higher[lower]) == 2
            theta *= Fraction(1, 2)
        # Every other nonroot has no higher neighbour across the cut.
        for v in range(n * n):
            if v not in {y * n + k - 1 for y in range(n)} | {y * n + k for y in range(n)}:
                assert all((u % n < k) == (v % n < k) for u in higher[v])
        assert theta == Fraction(1, 2 ** n)
        checked += 1
    assert checked == comb(16, 8)
    return {"n": n, "k": k, "shuffles_checked": checked, "theta_each": "1/16"}


def main():
    ans = {
        "n4_exhaustive": exhaustive_n4(),
        "n5_n6_structural": [structural_check(n, k) for n in (5, 6) for k in range(2, n - 1)],
        "status": "passed",
    }
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("interface-theta-family-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
