#!/usr/bin/env python3
"""Exact A2/B_top2 scan for finite connected simple graphs.

Words are vertex orders.  A birth occurs when the newly appended vertex has
no neighbour among earlier vertices.  A2 counts words with exactly two births.
B_top2 counts words whose first two vertices are non-adjacent and for which
every later vertex is on the boundary of the already exposed set.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


def masks(n: int, edges: list[tuple[int, int]]) -> list[int]:
    adj = [0] * n
    for u, v in edges:
        if not (0 <= u < n and 0 <= v < n and u != v):
            raise ValueError((n, u, v))
        if u > v:
            u, v = v, u
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def is_connected(adj: list[int]) -> bool:
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        x = frontier
        while x:
            bit = x & -x
            v = bit.bit_length() - 1
            nxt |= adj[v]
            x -= bit
        nxt &= ~seen
        seen |= nxt
        frontier = nxt
    return seen == (1 << len(adj)) - 1


def exact_counts(adj: list[int]) -> tuple[int, int]:
    """Return (A2, B_top2), with arbitrary precision integer arithmetic."""
    n = len(adj)
    full = (1 << n) - 1
    size = full + 1

    # dp0/1/2[S]: number of orders of S with the stated number of births.
    dp0 = [0] * size
    dp1 = [0] * size
    dp2 = [0] * size
    dp0[0] = 1
    for s in range(size):
        a0, a1, a2 = dp0[s], dp1[s], dp2[s]
        if not (a0 or a1 or a2):
            continue
        missing = full ^ s
        x = missing
        while x:
            bit = x & -x
            v = bit.bit_length() - 1
            t = s | bit
            born = (adj[v] & s) == 0
            if born:
                dp1[t] += a0
                dp2[t] += a1
            else:
                dp0[t] += a0
                dp1[t] += a1
                dp2[t] += a2
            x -= bit
    a2_total = dp2[full]

    # h[S]: boundary-growth completions after the initial exposed set S.
    h = [0] * size
    h[full] = 1
    for s in range(full - 1, 0, -1):
        neigh = 0
        x = s
        while x:
            bit = x & -x
            neigh |= adj[bit.bit_length() - 1]
            x -= bit
        boundary = neigh & (full ^ s)
        total = 0
        x = boundary
        while x:
            bit = x & -x
            total += h[s | bit]
            x -= bit
        h[s] = total
    btop2 = 0
    for u in range(n):
        for v in range(n):
            if u != v and not (adj[u] >> v) & 1:
                btop2 += h[(1 << u) | (1 << v)]
    return a2_total, btop2


def path(n: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(n - 1)]


def cycle(n: int) -> list[tuple[int, int]]:
    return path(n) + [(0, n - 1)]


def ladder(rungs: int) -> list[tuple[int, int]]:
    edges = []
    for i in range(rungs - 1):
        edges += [(i, i + 1), (rungs + i, rungs + i + 1)]
    edges += [(i, rungs + i) for i in range(rungs)]
    return edges


def binary_tree(height: int) -> list[tuple[int, int]]:
    n = (1 << (height + 1)) - 1
    return [(i, 2 * i + j) for i in range((n - 1) // 2) for j in (1, 2)]


def regular_rooted_tree(max_degree: int, height: int) -> list[tuple[int, int]]:
    """Root degree max_degree; other internal vertices have max_degree-1 children."""
    if max_degree < 2 or height < 1:
        raise ValueError
    levels = [[0]]
    edges = []
    next_vertex = 1
    for level in range(height):
        children_per_parent = max_degree if level == 0 else max_degree - 1
        nxt = []
        for parent in levels[-1]:
            for _ in range(children_per_parent):
                edges.append((parent, next_vertex))
                nxt.append(next_vertex)
                next_vertex += 1
        levels.append(nxt)
    return edges


def comb(teeth: int) -> list[tuple[int, int]]:
    # A path spine, one pendant leaf at every spine vertex: maximum degree 3.
    return path(teeth) + [(i, teeth + i) for i in range(teeth)]


def barbell(cycle_size: int, bridge_edges: int) -> list[tuple[int, int]]:
    # Two cycles joined by a bridge path.  Each attachment point has degree 3.
    if cycle_size < 3 or bridge_edges < 1:
        raise ValueError
    edges = cycle(cycle_size)
    offset = cycle_size
    edges += [(offset + i, offset + (i + 1) % cycle_size) for i in range(cycle_size)]
    # bridge_edges edges from 0 to offset, with bridge_edges-1 new internal vertices
    last = 0
    for j in range(bridge_edges - 1):
        z = 2 * cycle_size + j
        edges.append((last, z) if last < z else (z, last))
        last = z
    edges.append((last, offset) if last < offset else (offset, last))
    return sorted(set(edges))


def lollipop(cycle_size: int, tail_edges: int) -> list[tuple[int, int]]:
    if cycle_size < 3 or tail_edges < 1:
        raise ValueError
    edges = cycle(cycle_size)
    last = 0
    for j in range(tail_edges):
        z = cycle_size + j
        edges.append((last, z) if last < z else (z, last))
        last = z
    return sorted(set(edges))


@dataclass(frozen=True)
class Family:
    name: str
    param: int
    edges: list[tuple[int, int]]


def families(limit: int) -> list[Family]:
    out: list[Family] = []
    for n in range(2, limit + 1):
        out.append(Family("path", n, path(n)))
    for n in range(3, limit + 1):
        out.append(Family("cycle", n, cycle(n)))
    for r in range(2, limit // 2 + 1):
        out.append(Family("ladder", r, ladder(r)))
        out.append(Family("comb", r, comb(r)))
    for h in range(1, 8):
        e = binary_tree(h)
        n = 1 + max(max(x) for x in e)
        if n <= limit:
            out.append(Family("binary_tree", h, e))
    for degree in (4, 5):
        for h in range(1, 5):
            e = regular_rooted_tree(degree, h)
            n = 1 + max(max(x) for x in e)
            if n <= limit:
                out.append(Family(f"regular_tree_D{degree}", h, e))
    for c in range(3, 8):
        for b in range(1, 7):
            e = barbell(c, b)
            n = 1 + max(max(x) for x in e)
            if n <= limit:
                out.append(Family(f"barbell_c{c}", b, e))
    for c in range(3, 10):
        for t in range(1, 9):
            e = lollipop(c, t)
            n = 1 + max(max(x) for x in e)
            if n <= limit:
                out.append(Family(f"lollipop_c{c}", t, e))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=18)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for fam in families(args.limit):
        n = 1 + max(max(e) for e in fam.edges)
        adj = masks(n, fam.edges)
        assert is_connected(adj)
        a2, b2 = exact_counts(adj)
        delta = max(x.bit_count() for x in adj)
        rows.append({
            "family": fam.name, "parameter": fam.param, "n": n,
            "edges": fam.edges, "max_degree": delta,
            "A2": str(a2), "Btop2": str(b2),
            "ratio_A2_over_Btop2": f"{a2}/{b2}",
            "ratio_over_max_degree": f"{a2}/{b2 * delta}",
        })
        print(f"{fam.name:14} p={fam.param:2} n={n:2} D={delta} A2/B={a2}/{b2}")
    payload = {"method": "exact subset DP", "rows": rows}
    args.output.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
