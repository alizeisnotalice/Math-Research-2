#!/usr/bin/env python3
"""Independent layer census for the n=9 completion instance.

This deliberately does not read the C++ mask layers.  Raw independent-set
counts use a seven-row transfer DP.  Diagonal-reflection fixed counts use the
independent-set polynomial of the reflection quotient graph, where each orbit
has weight one or two according to its original cardinality.
"""
from collections import defaultdict
from functools import cache
import json
from pathlib import Path


Q = 7
FORBIDDEN = {(0, 0)}


def independent_rows():
    rows = [x for x in range(1 << Q) if not (x & (x << 1))]
    dp = {0: [1]}
    for y in range(Q):
        nxt = defaultdict(lambda: [0] * (Q * Q + 1))
        for prev, poly in dp.items():
            for row in rows:
                if row & prev or (y == 0 and row & 1):
                    continue
                shift = row.bit_count()
                target = nxt[row]
                for k, value in enumerate(poly):
                    if value:
                        target[k + shift] += value
        dp = nxt
    ans = [0] * (Q * Q + 1)
    for poly in dp.values():
        for k, value in enumerate(poly):
            ans[k] += value
    return ans


def adjacent(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


def fixed_polynomial():
    # A quotient vertex is an orbit under (x,y) -> (y,x).  Any orbit with an
    # internal grid edge is forbidden from an invariant independent set.
    seen = set()
    orbits = []
    for y in range(Q):
        for x in range(Q):
            cell = (x, y)
            if cell in seen or cell in FORBIDDEN:
                continue
            orbit = {cell, (y, x)}
            seen |= orbit
            if orbit & FORBIDDEN:
                continue
            if any(adjacent(a, b) for a in orbit for b in orbit if a != b):
                continue
            orbits.append(tuple(sorted(orbit)))
    n = len(orbits)
    neigh = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if any(adjacent(a, b) for a in orbits[i] for b in orbits[j]):
                neigh[i] |= 1 << j
                neigh[j] |= 1 << i
    weights = [len(x) for x in orbits]

    @cache
    def rec(mask):
        if not mask:
            return (1,)
        v = max((i for i in range(n) if mask >> i & 1),
                key=lambda i: (neigh[i] & mask).bit_count())
        no = rec(mask & ~(1 << v))
        yes = rec(mask & ~(neigh[v] | (1 << v)))
        out = list(no) + [0] * weights[v]
        for k, value in enumerate(yes):
            out[k + weights[v]] += value
        return tuple(out)

    return list(rec((1 << n) - 1)), n


def main():
    raw = independent_rows()
    fixed, quotient_vertices = fixed_polynomial()
    fixed += [0] * (len(raw) - len(fixed))
    reps = [(a + b) // 2 for a, b in zip(raw, fixed)]
    assert all((a + b) % 2 == 0 for a, b in zip(raw, fixed))
    result = {
        "n": 9,
        "U_geometry": "7x7 grid with (0,0) removed",
        "method": "row-transfer raw counts plus reflection quotient recurrence",
        "quotient_vertices": quotient_vertices,
        "raw": raw,
        "fixed": fixed,
        "reps": reps,
        "raw_total": sum(raw),
        "fixed_total": sum(fixed),
        "rep_total": sum(reps),
        "max_reps": max(reps),
        "max_k": reps.index(max(reps)),
    }
    # These expected integer values are the audited n=9 target, not input
    # from the C++ layer files.
    assert result["raw_total"] == 877160008
    assert result["fixed_total"] == 185378
    assert result["rep_total"] == 438672693
    assert result["max_reps"] == 83899706 and result["max_k"] == 12
    path = Path(__file__).with_name("n9-reflection-layer-audit-python.json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
