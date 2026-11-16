#!/usr/bin/env python3
"""Small exact checks for the relative-order interface-gluing ledger.

The 2-by-4 rectangle is enumerated completely.  The path calculation checks
the closed fibre formula U(P_m,r)=binom(m-1,r) for m<=8.
"""

from collections import Counter
from itertools import permutations
from math import comb
import json
from pathlib import Path
import sys


def rectangle(width, height):
    adj = []
    for y in range(height):
        for x in range(width):
            ns = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < width and 0 <= yy < height:
                    ns.append(yy * width + xx)
            adj.append(tuple(ns))
    return adj


def unique_peak_counts(adj):
    counts = Counter()
    for lab in permutations(range(len(adj))):
        roots = [v for v, ns in enumerate(adj) if all(lab[v] > lab[u] for u in ns)]
        if len(roots) == 1:
            counts[roots[0]] += 1
    return counts


def main():
    path_checks = {}
    for vertices in range(2, 9):
        got = unique_peak_counts(rectangle(vertices, 1))
        expected = {r: comb(vertices - 1, r) for r in range(vertices)}
        assert dict(got) == expected
        path_checks[str(vertices)] = expected

    # A 2-by-4 side of a vertical cut in P_4 square P_4.  The safe root
    # locations are the outer (x=0) column; their regional counts are exact.
    strip = unique_peak_counts(rectangle(2, 4))
    outer = {str(2 * y): strip[2 * y] for y in range(4)}
    assert sum(strip.values()) == 3584
    ans = {
        "path_unique_peak_formula_checked_through_vertices": 8,
        "path_counts": path_checks,
        "two_by_four_unique_peak_total": sum(strip.values()),
        "two_by_four_outer_column_safe_root_counts": outer,
        "status": "passed",
    }
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("interface-gluing-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
