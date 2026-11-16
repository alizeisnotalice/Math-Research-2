#!/usr/bin/env python3
"""Exact binary classification for H_d(x^2)H_d(x^3), d>=6.

For S={0,1,2,d-1,d}, collisions are equivalent to 3t and 2t both
belonging to the positive difference set of S.  The script checks the
closed-form classification against direct coordinate enumeration.
"""


def support(d):
    return {0, 1, 2, d - 1, d}


def directly_binary(d):
    s = support(d)
    values = [2*i + 3*j for i in s for j in s]
    return len(set(values)) == 25


def formula_binary(d):
    return d not in {6, 7, 9}


for d in range(6, 301):
    assert directly_binary(d) == formula_binary(d), d

print("H_d(x^2)H_d(x^3) is binary exactly for d>=6 with d not in {6,7,9}")
