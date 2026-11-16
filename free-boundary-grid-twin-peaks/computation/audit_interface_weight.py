#!/usr/bin/env python3
"""Exact n=2,n=3 audit of the natural-forest interface weight identity."""

from itertools import permutations
from fractions import Fraction
import json
from pathlib import Path
import sys


def graph(n):
    adj = []
    for y in range(n):
        for x in range(n):
            z = 0
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    z |= 1 << (yy * n + xx)
            adj.append(z)
    edges = [(v, u) for v, ns in enumerate(adj) for u in range(v + 1, n * n) if ns & (1 << u)]
    return adj, edges


def run(n):
    adj, edges = graph(n)
    N = n * n
    full = (1 << N) - 1
    checked_labels = checked_partitions = 0
    min_num, min_den = None, None
    max_cut = 0
    for lab in permutations(range(N)):
        hi = [sum((1 << u) for u in range(N) if adj[v] & (1 << u) and lab[u] > lab[v]) for v in range(N)]
        roots = [v for v in range(N) if hi[v] == 0]
        if len(roots) != 2:
            continue
        checked_labels += 1
        r0, r1 = roots
        free = [v for v in range(N) if v not in roots]
        for code in range(1 << len(free)):
            A = 1 << r0
            for j, v in enumerate(free):
                if code & (1 << j):
                    A |= 1 << v
            B = full ^ A
            assert B & (1 << r1)
            prod_a = prod_h = 1
            valid = True
            for v in range(N):
                if v == r0 or v == r1:
                    continue
                region = A if A & (1 << v) else B
                h = hi[v].bit_count()
                a = (hi[v] & region).bit_count()
                if a == 0:
                    valid = False
                    break
                prod_a *= a
                prod_h *= h
            if not valid:
                continue
            cut = sum(bool(A & (1 << v)) != bool(A & (1 << u)) for v, u in edges)
            checked_partitions += 1
            max_cut = max(max_cut, cut)
            # theta=prod_a/prod_h, and theta >= 2^(-cut).
            num, den = prod_a * (1 << cut), prod_h
            assert num >= den
            if min_num is None or num * min_den < min_num * den:
                min_num, min_den = num, den
    assert checked_labels == (8 if n == 2 else 143808)
    return {
        "n": n,
        "K2_labellings": checked_labels,
        "eligible_ordered_root_partitions": checked_partitions,
        "maximum_cut_size_seen": max_cut,
        "minimum_theta_times_2_to_cut": str(Fraction(min_num, min_den)),
        "status": "all tested partitions satisfy theta >= 2^(-cut)",
    }


def main():
    ans = {"n2": run(2), "n3": run(3)}
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("interface-weight-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
