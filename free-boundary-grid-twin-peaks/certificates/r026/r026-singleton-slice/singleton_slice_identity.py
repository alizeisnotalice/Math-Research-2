#!/usr/bin/env python3
"""Exact n=3 audit of the singleton first-merge slice identities.

The all-n proof is in SINGLETON_SLICE_IDENTITY.md.  This enumerates all 9!
words to check independently the two classifications used there: exact-one-
singleton first merges, lower-singleton (beta-positive) words, and the
global-top-two slice.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent


def graph(n: int) -> list[set[int]]:
    out = [set() for _ in range(n * n)]
    for v in range(n * n):
        x, y = v % n, v // n
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < n and 0 <= yy < n:
                out[v].add(yy * n + xx)
    return out


def components(adj: list[set[int]], seen: set[int]) -> list[set[int]]:
    rest, ans = set(seen), []
    while rest:
        stack, item = [rest.pop()], set()
        while stack:
            v = stack.pop()
            item.add(v)
            hit = adj[v] & rest
            rest -= hit
            stack.extend(hit)
        ans.append(item)
    return ans


def inspect(adj: list[set[int]], word: tuple[int, ...]):
    seen: set[int] = set()
    births: list[int] = []
    first = None
    for r, v in enumerate(word, 1):
        if not (adj[v] & seen):
            births.append(v)
        cs = components(adj, seen) if seen else []
        hit = [c for c in cs if c & adj[v]]
        if first is None and len(hit) >= 2:
            first = sorted(hit, key=lambda c: min(c))
        seen.add(v)
    return births, first


def main() -> None:
    n, adj = 3, graph(3)
    all_k2 = exact_one = beta_plus = top2_exact_one = 0
    base_mass = Fraction()
    q_hist = defaultdict(int)
    for word in itertools.permutations(range(n * n)):
        births, cs = inspect(adj, word)
        if len(births) != 2 or cs is None or len(cs) != 2:
            continue
        all_k2 += 1
        sizes = tuple(sorted(map(len, cs)))
        if sizes[0] != 1 or sizes[1] == 1:
            continue
        q = sum(sizes)
        exact_one += 1
        q_hist[q] += 1
        base_mass += Fraction(1, q)
        low_root = births[1]
        if any(low_root in c and len(c) == 1 for c in cs):
            beta_plus += 1
        if word[0] in cs[0] and word[1] in cs[1] or word[0] in cs[1] and word[1] in cs[0]:
            top2_exact_one += 1
    assert all_k2 == 143808
    # Each fixed exact-one-singleton fibre has q interleavings, q-1 of
    # which make its singleton root the lower root, and exactly two of
    # which lie in the global-top-two slice.
    assert beta_plus == exact_one - base_mass
    assert top2_exact_one == 2 * base_mass
    out = {
        "kind": "exact-n3-singleton-first-merge-slice-audit",
        "all_K2_words": all_k2,
        "exactly_one_singleton_first_merge_words_U": exact_one,
        "base_fibre_mass_sum_P": str(base_mass),
        "lower_singleton_words_Bplus_exact_one": beta_plus,
        "global_top_two_words_B1_exact_one": top2_exact_one,
        "q_histogram_among_full_words": dict(sorted(q_hist.items())),
        "identities": {
            "U_equals_Bplus_plus_P": True,
            "B1_equals_2P": True,
            "U_equals_Bplus_plus_B1_over_2": True,
        },
        "scope": "Exact n=3 audit only. The equalities themselves follow for every finite connected graph after grouping complete first-merge fibres; they exclude the separate two-singleton q=2 fibres.",
    }
    text = json.dumps(out, indent=2) + "\n"
    HERE.joinpath("singleton_slice_identity_n3.json").write_text(text)
    print(text, end="")
    print("output_sha256=" + hashlib.sha256(text.encode()).hexdigest())


if __name__ == "__main__":
    main()
