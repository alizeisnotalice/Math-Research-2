#!/usr/bin/env python3
"""Exact finite index audit for the b=11 and b>=12 R_8 closures."""

from math import gcd

S = {0, 1, 2, 3, 4, 5, 8}


def pairs(x, y, n):
    return {(i, n-i) for i in x if n-i in y}


for b, avals in ((11, range(2, 9)), *[(b, tuple(a for a in range(2, 9) if gcd(a, b) == 1)) for b in range(12, 41)]):
    for a in avals:
        assert gcd(a, b) == 1
        target = {a*i+b*j for i in S for j in S}
        assert len(target) == 49
        actual = []
        for i in range(17):
            for j in range(17):
                hit = [(r, s) for r in S for s in S if a*i+b*j == a*r+b*s]
                if hit and hit[0] != (i, j):
                    actual.append(((i, j), hit[0]))
        expected = [((i+b, j-a), (i, j)) for i in S if i+b <= 16 for j in S if j >= a]
        assert actual == expected

# Rows used after the standard low-prefix support reductions.
assert pairs(S, {0, 8}, 9) == {(1, 8)}                 # I
assert pairs({0, 1, 8}, {0, 1, 2, 3, 4, 8}, 9) == {(1, 8), (8, 1)}
assert pairs({0, 1, 8}, {0, 2, 3, 4}, 10) == {(8, 2)} # III after K9
assert pairs({0, 1, 2, 5}, {0, 1, 2, 3, 8}, 9) == {(1, 8)}
assert pairs({0, 1, 2, 5}, {0, 1, 2, 3}, 8) == {(5, 3)} # II0 after K9
x, y = {0, 1, 2, 5, 8}, {0, 1, 2, 3}
assert pairs(x, y, 9) == {(8, 1)}
assert pairs(x, {0, 2, 3}, 10) == {(8, 2)}
assert pairs(x, {0, 3}, 11) == {(8, 3)}
assert all(not pairs(x, {0, 3}, n) for n in range(12, 17))
assert pairs(x, {0, 3}, 1) == {(1, 0)}
assert pairs(x, {0, 3}, 3) == {(0, 3)}
assert pairs(x, {0, 3}, 4) == {(1, 3)}
print("b=11 and b=12..40 index/row assertions passed")
