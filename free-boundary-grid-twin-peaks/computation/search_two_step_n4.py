#!/usr/bin/env python3
"""Exact n=4 search for the completion-kernel two-step probability p(S).

For a connected high prefix S, h(S) is the number of legal connected
completion orders.  This computes
 p(S)=sum_{x in dS, y in N(x)\closed(S)} h(S+x+y)/h(S)
and tests the pointwise candidate p(S)<=Delta/|dS|.  It is a finite
counterexample/candidate audit only, not an asymptotic probability estimate.
"""
from functools import lru_cache
import json
from fractions import Fraction
from pathlib import Path

n = 4
N = n * n
FULL = (1 << N) - 1
adj = [0] * N
for v in range(N):
    x, y = v % n, v // n
    for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
        if 0 <= X < n and 0 <= Y < n:
            adj[v] |= 1 << (Y * n + X)


def neigh(mask):
    ans = 0
    while mask:
        b = mask & -mask
        mask ^= b
        ans |= adj[b.bit_length() - 1]
    return ans


@lru_cache(maxsize=None)
def h(mask):
    if mask == FULL:
        return 1
    bd = neigh(mask) & ~mask
    total = 0
    while bd:
        b = bd & -bd
        bd ^= b
        total += h(mask | b)
    return total


def connected(mask):
    if not mask:
        return False
    seen = mask & -mask
    frontier = seen
    while frontier:
        q = frontier
        frontier = 0
        while q:
            b = q & -q
            q ^= b
            frontier |= adj[b.bit_length() - 1]
        frontier &= mask & ~seen
        seen |= frontier
    return seen == mask


def cells(mask):
    return [[v % n, v // n] for v in range(N) if mask >> v & 1]


def main():
    best = None
    violations = []
    tested = 0
    for S in range(1, FULL):
        if not connected(S):
            continue
        bd = neigh(S) & ~S
        B = bd.bit_count()
        if not B:
            continue
        tested += 1
        closed = S | neigh(S)
        numerator = 0
        q = bd
        while q:
            xb = q & -q
            q ^= xb
            x = xb.bit_length() - 1
            ys = adj[x] & ~closed
            while ys:
                yb = ys & -ys
                ys ^= yb
                numerator += h(S | xb | yb)
        den = h(S)
        f = Fraction(numerator, den)
        rec = (f, B, S, numerator, den)
        if best is None or f * B > best[0] * best[1]:
            best = rec
        if numerator * B > 4 * den:
            violations.append(rec)
    f, B, S, num, den = best
    out = {
        "n": n, "tested_connected_prefixes": tested,
        "max_p": {"numerator": f.numerator, "denominator": f.denominator,
                  "boundary": B, "S_cells": cells(S),
                  "p_times_boundary": {"numerator": (f * B).numerator,
                                       "denominator": (f * B).denominator}},
        "candidate_p_le_4_over_boundary_violations": len(violations),
        "first_violation": (None if not violations else {
            "p": [violations[0][0].numerator, violations[0][0].denominator],
            "boundary": violations[0][1], "S_cells": cells(violations[0][2])}),
    }
    target = Path(__file__).parent / "results/two-step-n4-r006.json"
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
