#!/usr/bin/env python3
"""Exact n=3 check of T_{N-3}=8 A1(G-Nbar(corner))."""
from __future__ import annotations

import json
from itertools import permutations
from pathlib import Path


def one_peak_count(vertices, n=3):
    V = list(vertices)
    S = set(V)
    adj = {a: set() for a in V}
    for a in V:
        x, y = a % n, a // n
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            b = yy * n + xx
            if 0 <= xx < n and 0 <= yy < n and b in S:
                adj[a].add(b)
    ans = 0
    for word in permutations(V):
        pos = {a: i for i, a in enumerate(word)}
        ans += sum(all(pos[a] < pos[b] for b in adj[a]) for a in V) == 1
    return ans


def main(tail: Path, output: Path):
    x = json.loads(tail.read_text())
    n, N = x["n"], x["n"] ** 2
    assert n == 3
    global_hist = [sum(z["histogram_L"][i] for z in x["vertices"])
                   for i in range(N - 1)]
    a1 = one_peak_count(set(range(N)) - {0, 1, n})
    assert global_hist[N - 4] == 8 * a1
    output.write_text(json.dumps({"status": "passed", "n": n,
                                  "terminal_L": N - 3,
                                  "global_terminal_tail": global_hist[N - 4],
                                  "A1_corner_puncture": a1,
                                  "eight_times_A1": 8 * a1}, indent=2) + "\n")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(); p.add_argument("tail", type=Path); p.add_argument("output", type=Path)
    a = p.parse_args(); main(a.tail, a.output)
