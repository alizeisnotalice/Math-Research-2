#!/usr/bin/env python3
"""Independent integer interleaving DP for peak counts on trees.

Ranks are high-to-low.  The A2 recursion does not use the integral-message
implementation used by the comparison calculation.  For Btop2 we fix the
ordered two highest peak roots, delete them, and count no-peak relative
orders in every marked forest component by a separate integer recursion.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from math import comb, factorial


def full_binary(h: int) -> list[set[int]]:
    n = (1 << (h + 1)) - 1
    a = [set() for _ in range(n)]
    for v in range(1, n):
        p = (v - 1) // 2
        a[v].add(p)
        a[p].add(v)
    return a


def joined_binary(h: int, k: int) -> list[set[int]]:
    left = full_binary(h)
    right = full_binary(k)
    off = len(left)
    g = [set(x) for x in left] + [{off + w for w in x} for x in right]
    g[0].add(off)
    g[off].add(0)
    return g


def rooted_children(adj: list[set[int]], root: int, vertices: set[int] | None = None):
    if vertices is None:
        vertices = set(range(len(adj)))
    parent = {root: -1}
    order = [root]
    for v in order:
        for w in adj[v]:
            if w in vertices and w != parent[v]:
                parent[w] = v
                order.append(w)
    assert len(parent) == len(vertices)
    child = {v: [] for v in vertices}
    for v, p in parent.items():
        if p >= 0:
            child[p].append(v)
    return child


def merge_states(left, right, peak_limit: int):
    """Merge root-containing left block with a child subtree right block.

    A state (rank, k, h) records the root's rank in its block, peaks other
    than the root, and whether some child root precedes it.  The child root
    becomes a peak iff it precedes this root and had no earlier child.
    """
    out = defaultdict(int)
    for (m, r, k, h), av in left.items():
        for (c, q, kc, hc), bv in right.items():
            for s in range(c + 1):
                g = q <= s
                kp = k + kc + int(g and not hc)
                if kp > peak_limit:
                    continue
                ways = comb(r - 1 + s, s) * comb(m - r + c - s, c - s)
                out[(m + c, r + s, kp, h or g)] += av * bv * ways
    return out


def all_peak_states(adj: list[set[int]], root: int = 0, peak_limit: int = 2):
    child = rooted_children(adj, root)

    def rec(v):
        state = {(1, 1, 0, False): 1}
        for w in child[v]:
            state = merge_states(state, rec(w), peak_limit)
        return state

    return rec(root)


def count_a1_a2(adj: list[set[int]]):
    st = all_peak_states(adj)
    a1 = a2 = 0
    for (_, _, k, higher), ways in st.items():
        peaks = k + int(not higher)
        if peaks == 1:
            a1 += ways
        if peaks == 2:
            a2 += ways
    return a1, a2


def nopeak_component(adj, vertices: set[int], marked: set[int]):
    """Integer count of relative orders with no unmarked peak.

    A marked vertex has an already-higher external neighbour and therefore
    cannot be a peak.  Root at one marked vertex; at every merge an unmarked
    child root would be a peak exactly in the same g and hc situation as the
    A2 recursion, and is rejected.
    """
    root = min(marked)
    child = rooted_children(adj, root, vertices)

    def rec(v):
        # (size, rank, has-higher-child) -> integer count
        state = {(1, 1, False): 1}
        for w in child[v]:
            right = rec(w)
            out = defaultdict(int)
            for (m, r, h), av in state.items():
                for (c, q, hc), bv in right.items():
                    for s in range(c + 1):
                        g = q <= s
                        if w not in marked and g and not hc:
                            continue
                        ways = comb(r - 1 + s, s) * comb(m - r + c - s, c - s)
                        out[(m + c, r + s, h or g)] += av * bv * ways
            state = out
        return state

    return sum(rec(root).values())


def forest_after(adj: list[set[int]], removed: set[int]):
    rest = set(range(len(adj))) - removed
    out = []
    while rest:
        q = [next(iter(rest))]
        comp = {q[0]}
        rest.remove(q[0])
        for v in q:
            for w in adj[v]:
                if w in rest:
                    rest.remove(w)
                    comp.add(w)
                    q.append(w)
        out.append(comp)
    return out


def count_btop2(adj: list[set[int]]):
    n = len(adj)
    total = 0
    cache = {}
    for x in range(n):
        for y in range(n):
            if x == y or y in adj[x]:
                continue
            product = 1
            sizes = []
            for comp in forest_after(adj, {x, y}):
                marked = {v for v in comp if x in adj[v] or y in adj[v]}
                assert marked
                # Using explicit original vertex sets makes this cache merely
                # an optimisation within this independent integer program.
                key = (tuple(sorted(comp)), tuple(sorted(marked)))
                if key not in cache:
                    cache[key] = nopeak_component(adj, comp, marked)
                product *= cache[key]
                sizes.append(len(comp))
            interleave = factorial(n - 2)
            for size in sizes:
                interleave //= factorial(size)
            total += product * interleave
    return total


def hook_a1(adj: list[set[int]]):
    n = len(adj)
    ans = 0
    for root in range(n):
        child = rooted_children(adj, root)
        sizes = {}
        def size(v):
            sizes[v] = 1 + sum(size(w) for w in child[v])
            return sizes[v]
        size(root)
        den = 1
        for v in range(n):
            if v != root:
                den *= sizes[v]
        assert factorial(n - 1) % den == 0
        ans += factorial(n - 1) // den
    return ans
