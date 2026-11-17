#!/usr/bin/env python3
"""Exact n=3 audit of the two-row strip injection.

This enumerates the six R-completions, implements the proposed injection
into T-completions, and checks its cardinality and injectivity directly.
All counts are integers; no random sampling is used.
"""
from __future__ import annotations

import itertools
import json
from math import factorial
from pathlib import Path

n = 3
N = n * n
T = tuple(range(n))
B = tuple(range(n, 2 * n))
H = tuple(range(2 * n, N))
R = set(T + B)
Tset = set(T)


def nbrs(v: int) -> set[int]:
    x, y = v % n, v // n
    out: set[int] = set()
    for xx, yy in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
        if 0 <= xx < n and 0 <= yy < n:
            out.add(yy * n + xx)
    return out


def is_completion(initial: set[int], word: tuple[int, ...]) -> bool:
    seen = set(initial)
    for v in word:
        if v in seen or not (nbrs(v) & seen):
            return False
        seen.add(v)
    return True


def restricted_roots(sigma: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    roots: list[int] = []
    for v in sigma:
        if not (nbrs(v) & seen):
            roots.append(v)
        seen.add(v)
    return tuple(roots)


def shuffles_preserving_sigma(sigma: tuple[int, ...], free: tuple[int, ...]):
    """All words whose sigma-subsequence is sigma and whose other letters are free."""
    for word in itertools.permutations(sigma + free):
        if tuple(v for v in word if v in set(sigma)) == sigma:
            yield word


sigmas = [s for s in itertools.permutations(H) if is_completion(R, s)]
all_T = [w for w in itertools.permutations(tuple(v for v in range(N) if v not in Tset))
         if is_completion(Tset, w)]
constructed: dict[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = {}
per_sigma = []
for sigma in sigmas:
    P = restricted_roots(sigma)
    assert all(v // n == 2 for v in P)
    assert all(not (nbrs(x) & set(P)) for x in P)
    gates = tuple(v - n for v in P)
    free = tuple(v for v in B if v not in gates)
    count = 0
    for gate_order in itertools.permutations(gates):
        for tail in shuffles_preserving_sigma(sigma, free):
            word = gate_order + tail
            assert is_completion(Tset, word)
            assert word not in constructed, (word, constructed[word], sigma)
            constructed[word] = (sigma, gate_order, tail)
            count += 1
    M, p, q = len(H), len(P), len(free)
    expected = factorial(p) * factorial(M + q) // factorial(M)
    assert count == expected
    per_sigma.append({"sigma": list(sigma), "restricted_roots": list(P),
                      "p": p, "constructed": count, "formula": expected})

assert len(sigmas) == 6
assert len(all_T) == 156
assert len(constructed) == sum(x["constructed"] for x in per_sigma) == 96
assert set(constructed).issubset(set(all_T))
out = {
    "scope": "Exact n=3 audit of the strip injection; not an asymptotic proof by computation.",
    "n": n,
    "H_R": len(sigmas),
    "H_T": len(all_T),
    "injection_images": len(constructed),
    "per_sigma": per_sigma,
    "bound_using_maximum_p": 6 * 8,
    "formula_minimum_C_p": 8,
    "all_constructed_words_distinct": True,
}
Path(__file__).with_name("strip_entropy_n3.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
