#!/usr/bin/env python3
"""Exact audit for the L1-layer leaf-splitting parent tree on P_n square P_n.

At L1 layer k=x+y, use a_k=2(k-2^floor(log2 k)).  The parent map splits
exactly leaf a_k of the preceding breadth-first, left-to-right binary leaf
tree.  The finite square is a restriction of the first-quadrant tree.
All checks below use integer arithmetic.
"""

from collections import defaultdict


def split_index(k):
    assert k >= 1
    return 2 * (k - (1 << (k.bit_length() - 1)))


def parent(v):
    x, y = v
    if (x, y) == (0, 0):
        return None
    k = x + y
    a = split_index(k)
    if x <= a:
        return (x, y - 1)
    return (x - 1, y)


def audit(n):
    vertices = [(x, y) for x in range(n) for y in range(n)]
    vertex_set = set(vertices)
    children = defaultdict(list)
    for v in vertices:
        if v == (0, 0):
            continue
        p = parent(v)
        assert p in vertex_set
        assert p[0] + p[1] == v[0] + v[1] - 1
        children[p].append(v)

    # Strict decrease of x+y proves acyclicity and every vertex reaches root.
    subtree = {v: 1 for v in vertices}
    for v in sorted(vertices, key=lambda z: z[0] + z[1], reverse=True):
        p = parent(v)
        if p is not None:
            subtree[p] += subtree[v]

    maximum_by_layer = {}
    for r in range(2 * n - 1):
        layer = [v for v in vertices if v[0] + v[1] == r]
        if not layer:
            continue
        # Infinite-quadrant leaf bound, summed over the finite square's layers.
        layer_bound_numer = 4 * sum(k + 1 for k in range(r, 2 * n - 1))
        for v in layer:
            assert subtree[v] * (r + 1) <= layer_bound_numer, (n, v, subtree[v])
        maximum_by_layer[r] = max(subtree[v] for v in layer)

    return n, maximum_by_layer


def audit_leaf_schedule(limit):
    """Directly check the breadth-first leaf-splitting invariant."""
    leaves = [0]  # leaf depths at layer 0, in left-to-right order
    for k in range(1, limit + 1):
        a = split_index(k)
        assert 0 <= a < len(leaves), (k, a, len(leaves))
        depth = leaves[a]
        leaves[a:a + 1] = [depth + 1, depth + 1]
        lower = (k + 1).bit_length() - 1
        assert min(leaves) >= lower, (k, leaves)
        assert max(leaves) <= lower + 1, (k, leaves)


if __name__ == "__main__":
    audit_leaf_schedule(8192)
    for n in (8, 16, 32, 64, 128):
        n, mx = audit(n)
        r = n // 2
        # For r>=n/2, the theorem's coarse 16n bound follows from 8n^2/(r+1).
        assert max(value for level, value in mx.items() if level >= r) <= 16 * n
        print({"n": n, "max_at_or_beyond_n_over_2":
               max(value for level, value in mx.items() if level >= r),
               "at_layer": r, "layer_max": mx[r]})
