#!/usr/bin/env python3
"""Finite symbolic-index audit of the H_d carry case split (not a proof by samples)."""

from math import gcd

def active_rows(d, a):
    S = {0, 1, 2, d-1, d}
    positive = {j-a for j in S if j >= a}
    negative = {j for j in S if j >= a}
    return positive, negative

for d in range(6, 61):
    for a in range(2, d+1):
        pos, neg = active_rows(d, a)
        if a == 2 and d >= 7:
            assert {d-3, d-2} <= pos
            assert {d-3, d-2}.isdisjoint({0, 1, 2, 3, d-1, d})
        if a >= 3:
            s = d-1-a
            if s >= 4:
                assert {s, s+1} <= pos
                assert {s, s+1}.isdisjoint({0, 1, 2, 3, d-1, d})
            if s == 3:
                assert pos == {3, 4}
            if s == 2:
                assert pos == {2, 3}
            if s == 1:
                assert pos == {1, 2}
            if s == 0:
                assert pos == {0, 1}
            if s == -1:
                assert pos == {0}
print("symbolic-index case assertions passed for d=6..60")

# Independent finite check of the congruence bookkeeping used to derive (N).
# The written proof, not this sample, handles arbitrary d,a,b.
for d in range(6, 31):
    S = {0, 1, 2, d - 1, d}
    for a in range(1, d + 1):
        for b in range(d + 1, 2 * d + 1):
            if gcd(a, b) != 1:
                continue
            targets = {i * a + j * b: (i, j) for i in S for j in S}
            assert len(targets) == len(S) ** 2
            for i in range(2 * d + 1):
                for j in range(2 * d + 1):
                    exponent = i * a + j * b
                    if exponent not in targets:
                        continue
                    i0, j0 = targets[exponent]
                    assert (i - i0) % b == 0
                    t = (i - i0) // b
                    assert j == j0 - t * a
                    assert t in {0, 1}
                    if t == 1:
                        assert i0 <= 2 * d - b < d
                        assert j0 >= a
print("finite lift/carry bookkeeping checks passed for d=6..30, d<b<=2d")
