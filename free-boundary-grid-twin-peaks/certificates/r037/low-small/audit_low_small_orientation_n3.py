#!/usr/bin/env python3
"""Exact P3xP3 audit of the low-/high-small branch orientation identity."""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations
import hashlib
import json
import platform
from pathlib import Path

HERE = Path(__file__).resolve().parent


def grid(n=3):
    adj = [set() for _ in range(n * n)]
    for y in range(n):
        for x in range(n):
            u = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    adj[u].add(yy * n + xx)
    return adj


def components(seen, adj):
    left = set(seen)
    ans = []
    while left:
        c = {left.pop()}
        todo = list(c)
        while todo:
            u = todo.pop()
            hit = adj[u] & left
            left -= hit
            c |= hit
            todo.extend(hit)
        ans.append(c)
    return ans


def births_and_first_merge(word, adj):
    seen = set()
    births = []
    first_hit = None
    for v in word:
        cs = components(seen, adj) if seen else []
        hit = [c for c in cs if c & adj[v]]
        if not hit:
            births.append(v)
        if first_hit is None and len(hit) >= 2:
            first_hit = hit
        seen.add(v)
    return births, first_hit


def geometry(u, v, n=3):
    dx = abs(u % n - v % n)
    dy = abs(u // n - v // n)
    return dx + dy, Fraction(dx * dx + dy * dy, n * n)


def pair_to_json(x):
    assert isinstance(x, Fraction)
    return [x.numerator, x.denominator]


def main():
    adj = grid()
    low = defaultdict(lambda: [Fraction(), Fraction(), Fraction()])
    high_weighted = defaultdict(lambda: [Fraction(), Fraction(), Fraction()])
    total_k2 = 0
    for word in permutations(range(9)):
        births, hit = births_and_first_merge(word, adj)
        if len(births) != 2:
            continue
        assert hit is not None
        assert len(hit) == 2
        total_k2 += 1
        hi, lo = births
        high_branch = next(c for c in hit if hi in c)
        low_branch = next(c for c in hit if lo in c)
        a, b = len(high_branch), len(low_branch)
        q = a + b
        d, e = geometry(hi, lo)
        fs = (Fraction(1), Fraction(d), e)
        for j, f in enumerate(fs):
            low[b][j] += f
            high_weighted[a][j] += Fraction(b, a) * f

    assert sum(v[0] for v in low.values()) == total_k2
    assert set(low) == set(high_weighted)
    rows = {}
    for b in sorted(low):
        assert low[b] == high_weighted[b]
        rows[str(b)] = {
            "low_branch_raw": [pair_to_json(x) for x in low[b]],
            "high_branch_weighted_by_other_over_b": [
                pair_to_json(x) for x in high_weighted[b]
            ],
        }
    result = {
        "kind": "exact-P3xP3-low-small-orientation",
        "python": platform.python_version(),
        "K2_word_count": total_k2,
        "observables": ["1", "Manhattan distance", "squared Euclidean distance / 9"],
        "identity_checked": (
            "sum_{L_b} F = sum_{H_b} ((q-b)/b) F; "
            "here H_b has high branch size b"
        ),
        "rows_by_exact_small_branch_size": rows,
        "all_identities_passed": True,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    HERE.joinpath("audit-low-small-orientation-n3.json").write_text(text)
    print(text, end="")
    print("output_sha256=" + hashlib.sha256(text.encode()).hexdigest())


if __name__ == "__main__":
    main()
