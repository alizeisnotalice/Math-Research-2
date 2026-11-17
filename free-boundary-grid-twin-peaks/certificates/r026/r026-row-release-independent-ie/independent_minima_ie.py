#!/usr/bin/env python3
"""Independent exact check using inclusion-exclusion over local minima.

For a variable graph Q with an earlier fixed boundary, a legal reveal word is
exactly a permutation in which no forbidden vertex is a local minimum.  For an
independent event set S, requiring all v in S to be local minima imposes only
the height-two poset v < u for every Q-edge vu.  Its linear extensions are
counted by an order-ideal DP on S (never on all vertices of Q).
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter


@dataclass(frozen=True)
class Case:
    name: str
    n: int
    vertices: tuple[tuple[int, int], ...]
    forbidden: tuple[int, ...]
    adj: tuple[int, ...]


def make_case(n: int, kind: str) -> Case:
    # HT: variables B union H (all but initial top row); forbid minima in H.
    # HR: variables H; the initial B makes y=2 automatically non-minimal.
    if kind == "HT":
        vertices = tuple((x, y) for y in range(1, n) for x in range(n))
        forbid_coord = lambda x, y: y >= 2
    elif kind == "HR":
        vertices = tuple((x, y) for y in range(2, n) for x in range(n))
        forbid_coord = lambda x, y: y >= 3
    else:
        raise ValueError(kind)
    pos = {v: i for i, v in enumerate(vertices)}
    adj = []
    for x, y in vertices:
        mask = 0
        for xx, yy in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
            j = pos.get((xx, yy))
            if j is not None:
                mask |= 1 << j
        adj.append(mask)
    forbidden = tuple(i for i, (x, y) in enumerate(vertices) if forbid_coord(x, y))
    return Case(kind, n, vertices, forbidden, tuple(adj))


def independent_sets(case: Case):
    """Yield masks in the compact forbidden-coordinate index set, by DFS."""
    f = case.forbidden
    fadj = []
    for i, v in enumerate(f):
        z = 0
        for j, w in enumerate(f):
            if case.adj[v] >> w & 1:
                z |= 1 << j
        fadj.append(z)

    def rec(i: int, chosen: int, blocked: int):
        if i == len(f):
            yield chosen
            return
        yield from rec(i + 1, chosen, blocked)
        if not (blocked >> i & 1):
            yield from rec(i + 1, chosen | (1 << i), blocked | fadj[i] | (1 << i))
    yield from rec(0, 0, 0)


def le_height_two(case: Case, compact_S: int) -> int:
    """#linear extensions of the poset 'each selected v precedes its neighbors'."""
    selected = [case.forbidden[j] for j in range(len(case.forbidden)) if compact_S >> j & 1]
    m = len(selected)
    N = len(case.vertices)
    if not m:
        # Explicit factorial avoids calling the DP in the dominant empty case.
        z = 1
        for q in range(2, N + 1):
            z *= q
        return z
    source_bit_of_vertex = {v: 1 << j for j, v in enumerate(selected)}
    need = [0] * N
    for u in range(N):
        for v, b in source_bit_of_vertex.items():
            if case.adj[v] >> u & 1:
                need[u] |= b
    # Independence ensures selected sources themselves have no predecessors.
    assert all(need[v] == 0 for v in selected)
    nonsources = [u for u in range(N) if u not in source_bit_of_vertex]
    q = len(nonsources)
    e = [0] * (1 << m)
    for A in range(1 << m):
        e[A] = sum((need[u] & ~A) == 0 for u in nonsources)

    # dp[A][k]: exact number of prefix words containing exactly sources A and
    # k non-sources.  Existing non-sources remain eligible as A grows, hence
    # e[A]-k is exactly the number of possible next non-source vertices.
    dp = [[0] * (q + 1) for _ in range(1 << m)]
    dp[0][0] = 1
    full = (1 << m) - 1
    for A in range(1 << m):
        for k, value in enumerate(dp[A]):
            if not value:
                continue
            unused = full ^ A
            while unused:
                b = unused & -unused
                unused ^= b
                dp[A | b][k] += value
            if k < q:
                choices = e[A] - k
                if choices:
                    dp[A][k + 1] += value * choices
    return dp[full][q]


def estimate(case: Case) -> dict:
    hist: dict[int, int] = defaultdict(int)
    ideals = 0
    for S in independent_sets(case):
        m = S.bit_count()
        hist[m] += 1
        ideals += 1 << m
    N = len(case.vertices)
    # Largest live DP uses (2^m)*(N-m+1) Python integer slots.  This is a
    # conservative slot count, not a measured byte count.
    max_m = max(hist)
    peak_slots = (1 << max_m) * (N - max_m + 1)
    return {"case": case.name, "n": case.n, "variable_vertices": N,
            "forbidden_vertices": len(case.forbidden), "independent_set_count": sum(hist.values()),
            "independent_set_size_histogram": {str(k): hist[k] for k in sorted(hist)},
            "sum_order_ideal_states_upper": ideals,
            "max_selected_minima": max_m, "peak_dp_integer_slots_upper": peak_slots}


def compute(case: Case) -> dict:
    total = 0
    terms = 0
    pos = neg = 0
    t0 = perf_counter()
    for S in independent_sets(case):
        term = le_height_two(case, S)
        terms += 1
        if S.bit_count() & 1:
            total -= term
            neg += term
        else:
            total += term
            pos += term
    return {"case": case.name, "n": case.n, "value": str(total), "terms": terms,
            "positive_sum": str(pos), "negative_sum": str(neg),
            "elapsed_seconds": perf_counter() - t0}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--mode", choices=("estimate", "run"), required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for n in args.n:
        for kind in ("HT", "HR"):
            case = make_case(n, kind)
            rows.append(estimate(case) if args.mode == "estimate" else compute(case))
    args.output.write_text(json.dumps({"algorithm": "local-minimum-inclusion-exclusion-height-two-order-ideal", "mode": args.mode, "rows": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
