#!/usr/bin/env python3
"""Exact n=3 check of the single exceptional-birth fact used in r014."""

from itertools import permutations
import json
from pathlib import Path
import sys


def grid(n):
    out = []
    for y in range(n):
        for x in range(n):
            ns = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    ns.append(yy * n + xx)
            out.append(tuple(ns))
    return out


def main():
    n = 3
    adj = grid(n)
    k2 = 0
    birth_time_hist = {}
    for labels in permutations(range(n * n)):
        roots = [v for v in range(n * n) if all(labels[v] > labels[u] for u in adj[v])]
        if len(roots) != 2:
            continue
        order = sorted(range(n * n), key=labels.__getitem__, reverse=True)
        seen = {order[0]}
        exceptional = []
        for step, v in enumerate(order[1:], start=2):
            if not any(u in seen for u in adj[v]):
                exceptional.append(step)
            seen.add(v)
        assert len(exceptional) == 1
        assert order[exceptional[0] - 1] in roots
        assert order[0] in roots
        k2 += 1
        birth_time_hist[str(exceptional[0])] = birth_time_hist.get(str(exceptional[0]), 0) + 1
    assert k2 == 143808
    ans = {"n": n, "K2_labellings": k2, "every_K2_has_one_nonboundary_addition": True,
           "second_peak_birth_step_histogram": birth_time_hist, "status": "passed"}
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("two-peak-birth-audit-n3.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
