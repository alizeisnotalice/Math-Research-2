#!/usr/bin/env python3
"""Exact 3x3 check of the stable-frontload completion lower bound.

For every connected nonempty S, M(S) counts all one-peak orders whose root
lies in S, and h(S) counts legal completions after the already exposed set S.
The audited inequality is M(S) <= |S|(N-1)_(|S|-1) h(S).
"""
import json
import math
from pathlib import Path

N = 9
ALL = (1 << N) - 1
HERE = Path(__file__).resolve().parent


def is_connected(mask, adj):
    if not mask:
        return False
    seen = mask & -mask
    while True:
        boundary = 0
        bits = seen
        while bits:
            bit = bits & -bits
            v = bit.bit_length() - 1
            boundary |= adj[v]
            bits -= bit
        new = (boundary & mask) & ~seen
        if not new:
            return seen == mask
        seen |= new


def main():
    adj = [0] * N
    for v in range(N):
        x, y = v % 3, v // 3
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < 3 and 0 <= yy < 3:
                adj[v] |= 1 << (yy * 3 + xx)
    h = [0] * (1 << N)
    h[ALL] = 1
    for mask in range(ALL - 1, 0, -1):
        frontier = 0
        bits = mask
        while bits:
            bit = bits & -bits
            frontier |= adj[bit.bit_length() - 1]
            bits -= bit
        frontier &= ~mask
        while frontier:
            bit = frontier & -frontier
            h[mask] += h[mask | bit]
            frontier -= bit
    rooted = [h[1 << r] for r in range(N)]
    checks, minimum_slack = 0, None
    by_size = {}
    for mask in range(1, 1 << N):
        if not is_connected(mask, adj):
            continue
        d = mask.bit_count()
        m = sum(rooted[r] for r in range(N) if mask >> r & 1)
        rhs = d * math.perm(N - 1, d - 1) * h[mask]
        assert m <= rhs
        slack = rhs - m
        if minimum_slack is None or slack < minimum_slack:
            minimum_slack = slack
        entry = by_size.setdefault(d, {"sets": 0, "minimum_rhs_over_M": None})
        entry["sets"] += 1
        if m:
            ratio = rhs / m
            old = entry["minimum_rhs_over_M"]
            entry["minimum_rhs_over_M"] = ratio if old is None else min(old, ratio)
        checks += 1
    out = {"n": 3, "A1": str(sum(rooted)), "rooted_counts": rooted,
           "connected_sets_checked": checks, "minimum_integer_slack": str(minimum_slack),
           "by_size": by_size,
           "inequality": "M(S) <= |S|*(N-1)_(|S|-1)*h(S)"}
    (HERE / "frontload-n3-check.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
