#!/usr/bin/env python3
"""Independent exact recurrence on H=G-v for D4 representatives, n<=6.

This intentionally uses C_H(empty)=(N-1)! rather than the original full-grid
C_G table.  It is only a finite cross-check of the punctured formula.
"""
from __future__ import annotations

import gc
import argparse
import json
import sys
from collections import defaultdict
from math import factorial
from pathlib import Path


def graph(n: int):
    N = n * n
    adj = [0] * N
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    adj[v] |= 1 << (yy * n + xx)
    return adj, [a | (1 << v) for v, a in enumerate(adj)]


def independent_bags(n: int, allowed: int):
    """Independent sets contained in allowed, grouped by their size."""
    rows = [r for r in range(1 << n) if not r & (r << 1)]
    allowed_rows = [(allowed >> (n * y)) & ((1 << n) - 1) for y in range(n)]
    bags = defaultdict(list)

    def visit(y: int, previous: int, I: int):
        if y == n:
            bags[I.bit_count()].append(I)
            return
        for r in rows:
            if not (r & previous) and not (r & ~allowed_rows[y]):
                visit(y + 1, r, I | (r << (n * y)))

    visit(0, 0, 0)
    return bags


def exact_at_hole(n: int, v: int):
    """The stated H=G-v recurrence and its two alternating sums."""
    N = n * n
    adj, closed = graph(n)
    full = (1 << N) - 1
    Z = full & ~closed[v]
    bags = independent_bags(n, Z)
    C = {0: factorial(N - 1)}
    A = 0
    TB = C[0] * Z.bit_count()  # I=empty term in the B numerator.
    for k in range(1, max(bags) + 1):
        for I in bags[k]:
            cl = I
            q = I
            while q:
                bit = q & -q
                q ^= bit
                cl |= adj[bit.bit_length() - 1]
            numerator = 0
            q = I
            while q:
                bit = q & -q
                q ^= bit
                numerator += C[I ^ bit]
            den = cl.bit_count()
            assert numerator % den == 0, (n, v, I, numerator, den)
            count = numerator // den
            C[I] = count
            A += (1 if k & 1 else -1) * k * count
            avail = (full & ~(closed[v] | cl)).bit_count()
            TB += (-1 if k & 1 else 1) * count * avail
    assert TB % (N - 1) == 0, (n, v, TB)
    result = (A, TB // (N - 1), sum(map(len, bags.values())))
    del C, bags
    gc.collect()
    return result


def d4_orbit(n: int, x: int, y: int):
    m = n - 1
    out = set()
    for a, b in ((x, y), (y, x)):
        for sx in (a, m - a):
            for sy in (b, m - b):
                out.add((sx, sy))
    return out


def expected(n: int, r027: Path, r028: Path):
    if n <= 5:
        data = json.loads((r027 / f"punctured-table-n{n}.json").read_text())
        return [(x["A2_v"], x["B_v"]) for x in data["per_vertex"]]
    data = json.loads((r028 / "punctured-table-n6.json").read_text())
    return [(int(x["A2_v"]), int(x["B_v"])) for x in data["per_vertex"]]


def main(r027: Path, r028: Path):
    summary = []
    for n in range(3, 7):
        ref = expected(n, r027, r028)
        # Check the already stored table is D4 invariant before one
        # representative per orbit is used for the new recurrence.
        for y in range(n):
            for x in range(n):
                here = ref[y * n + x]
                assert all(ref[yy * n + xx] == here for xx, yy in d4_orbit(n, x, y))
        reps = [(x, y) for y in range((n + 1) // 2) for x in range(y + 1)]
        for x, y in reps:
            v = y * n + x
            got_a, got_b, sets = exact_at_hole(n, v)
            assert (got_a, got_b) == ref[v], (n, (x, y), (got_a, got_b), ref[v])
            summary.append({"n": n, "xy": [x, y], "independent_sets": sets,
                            "A2_v": got_a, "B_v": got_b})
            print(f"n={n}, v=({x},{y}), independent={sets}: matched", file=sys.stderr)
    print(json.dumps({"status": "passed", "representatives": summary}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--r027", type=Path, required=True,
                        help="Directory holding R027 punctured-table-n3..n5 JSON files")
    parser.add_argument("--r028", type=Path, required=True,
                        help="Directory holding R028 punctured-table-n6.json")
    args = parser.parse_args()
    main(args.r027, args.r028)
