#!/usr/bin/env python3
"""Exact root-resolved stable-fibre tails on punctured P_n square P_n.

For a deleted maximum v, and an exterior root p, this enumerates exactly the
words in D_{v,p}: their first vertex is p and every later vertex outside
N_G(v) has an earlier neighbour.  This is precisely the condition that p is
the unique peak outside the hole neighbourhood.  It returns the histogram of
L=J-1, where J is the first position of a neighbour of p.

The n=3,4 runs are deliberately a prefix DP on H=G-v, independent of the
independent-set inclusion-exclusion formula used for the stored A2_v,B_v
tables.  They are a finite check of the fibre identity, not an asymptotic
estimate.
"""
from __future__ import annotations

import argparse
import json
from functools import lru_cache
from pathlib import Path


def graph(n: int) -> list[int]:
    a = [0] * (n * n)
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    a[v] |= 1 << (yy * n + xx)
    return a


def tail_for_v(n: int, v: int) -> dict:
    """Return exact histogram and (B_v,A2_v) reconstructed from all p."""
    N = n * n
    full = (1 << N) - 1
    adj = graph(n)
    nv = adj[v]
    H = full ^ (1 << v)
    exterior = H & ~nv
    aggregate = [0] * N
    by_root = []

    q = exterior
    while q:
        bp = q & -q
        q ^= bp
        p = bp.bit_length() - 1

        @lru_cache(maxsize=None)
        def poly(S: int) -> tuple[int, ...]:
            """Polynomial in future non-neighbours before the first p-neighbour.

            The returned coefficient at r counts completions with r *additional*
            non-neighbour positions before the first p-neighbour.  Once a
            p-neighbour is in S the polynomial is the constant completion count.
            """
            if S == H:
                return (1,)
            p_seen_neighbour = bool(adj[p] & S)
            ans = [0] * (N - S.bit_count())
            missing = H & ~S
            r = missing
            while r:
                bx = r & -r
                r ^= bx
                x = bx.bit_length() - 1
                # A vertex outside N(v) must not be an additional peak.
                if not ((nv >> x) & 1) and not (adj[x] & S):
                    continue
                child = poly(S | bx)
                shift = (not p_seen_neighbour) and not ((adj[p] >> x) & 1)
                for i, c in enumerate(child):
                    ans[i + int(shift)] += c
            return tuple(ans)

        # p itself contributes one to L=J-1.
        hist0 = poly(bp)
        hist = [0] * N
        for r, c in enumerate(hist0):
            hist[r + 1] = c
            aggregate[r + 1] += c
        by_root.append({"p": p, "histogram_L": hist[1:],
                        "count": sum(hist), "sum_L": sum(i * c for i, c in enumerate(hist))})

    B = sum(aggregate)
    A = sum(i * c for i, c in enumerate(aggregate))
    return {"v": v, "xy": [v % n, v // n], "histogram_L": aggregate[1:],
            "B_v_from_D": B, "A2_v_from_fibres": A,
            "mean_L": {"numerator": A, "denominator": B}, "by_root": by_root}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int, choices=(3, 4))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = {"n": args.n, "method": "exact prefix DP on H=G-v", "vertices": []}
    for v in range(args.n * args.n):
        out["vertices"].append(tail_for_v(args.n, v))
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
