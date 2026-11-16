#!/usr/bin/env python3
"""Exact 3x3 audit of forced root-neighbour natural-forest edges."""

from collections import Counter
from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path
import sys


def grid(n):
    adj = []
    for y in range(n):
        for x in range(n):
            ns = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    ns.append(yy * n + xx)
            adj.append(tuple(ns))
    return adj


def main():
    n = 3
    adj = grid(n)
    far_count = 0
    total_force_mass = Fraction()
    min_force = None
    hist = Counter()
    distance_hist = Counter()
    max_forced_vertices = 0
    for lab in permutations(range(n * n)):
        higher = [tuple(u for u in adj[v] if lab[u] > lab[v]) for v in range(n * n)]
        roots = [v for v in range(n * n) if not higher[v]]
        if len(roots) != 2:
            continue
        p, q = roots
        dist = abs(p % n - q % n) + abs(p // n - q // n)
        if dist < 3:
            continue
        distance_hist[dist] += 1
        np, nq = set(adj[p]), set(adj[q])
        assert not (np & nq)
        forced = np | nq
        denominator = 1
        for v in forced:
            parent = p if v in np else q
            assert parent in higher[v]
            denominator *= len(higher[v])
        force_probability = Fraction(1, denominator)
        assert force_probability >= Fraction(1, 4 ** 8)
        far_count += 1
        total_force_mass += force_probability
        min_force = force_probability if min_force is None else min(min_force, force_probability)
        hist[str(force_probability)] += 1
        max_forced_vertices = max(max_forced_vertices, len(forced))
    # Independent complete distance ledger within the K=2 enumeration.
    assert distance_hist == Counter({3: 64512, 4: 15744})
    assert far_count == 80256
    assert total_force_mass >= Fraction(far_count, 4 ** 8)
    assert total_force_mass <= far_count
    ans = {
        "n": n,
        "K2_distance_at_least_3_labellings": far_count,
        "distance_histogram": dict(sorted(distance_hist.items())),
        "max_forced_root_neighbours": max_forced_vertices,
        "minimum_forced_event_probability": str(min_force),
        "sum_forced_event_probabilities": str(total_force_mass),
        "coarse_lower_bound": str(Fraction(far_count, 4 ** 8)),
        "probability_histogram": dict(sorted(hist.items(), key=lambda x: Fraction(x[0]))),
        "status": "passed",
    }
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("safe-root-forcing-n3.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
