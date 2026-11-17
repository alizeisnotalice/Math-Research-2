#!/usr/bin/env python3
"""Exact local-prefix witness for the degree-seven S={3} analysis.

This is deliberately not a global model: the displayed G,Q entries are only
low prefixes, and the script claims no all-root annulus or infinite extension.
It proves that rows through degree 7 plus Newton bounds through s_7 alone do
not yield a contradiction.
"""

import json
from fractions import Fraction
from pathlib import Path


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def scale(k, x): return (k * x[0], k * x[1])
def mul(x, y):
    # a=sqrt(2), so a^2=2.
    return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def conv(x, y, limit):
    out = [(0, 0)] * (limit + 1)
    for i, u in enumerate(x):
        for j, v in enumerate(y):
            if i + j <= limit:
                out[i + j] = add(out[i + j], mul(u, v))
    return out


def newton(P):
    # P=prod(1-w_i x); s_k+p1*s_{k-1}+...+p_{k-1}*s1+k*p_k=0.
    s = [(0, 0)]
    for k in range(1, len(P)):
        total = scale(-k, P[k])
        for j in range(1, k):
            total = add(total, scale(-1, mul(P[j], s[k-j])))
        s.append(total)
    return s


def main():
    z, one, a = (0, 0), (1, 0), (0, 1)
    P = [one, z, z, a, z, z, z, one]
    Pb = [one, z, z, scale(-1, a), z, z, z, one]
    # Formal prefixes only; they may be extended in higher degree.
    G = [one, z, z, (-1, 0), z, z, (2, 0), (-1, 0)]
    Q = [one, z, z, (2, 0), z, z, (3, 0), (-1, 0)]
    A, B = conv(P, G, 7), conv(Pb, Q, 7)
    C = conv(A, B, 7)
    assert A == [one, z, z, (-1, 1), z, z, (2, -1), z]
    assert B == [one, z, z, (2, -1), z, z, (3, -2), z]
    assert C == [one, z, z, one, z, z, one, z]
    s = newton(P)
    assert s == [z, z, z, (0, -3), z, z, (6, 0), (-7, 0)]
    # Exact strict checks: 3sqrt(2)<203/8, 6<12355/64, 7<49283/128.
    assert 18 * 64 < 203**2
    assert 6 * 64 < 12355 and 7 * 128 < 49283
    out = {
        "status": "passed_local_prefix_only",
        "P": [list(x) for x in P],
        "G_prefix_0_7": [list(x) for x in G],
        "Q_prefix_0_7": [list(x) for x in Q],
        "A_prefix_0_7": [list(x) for x in A],
        "B_prefix_0_7": [list(x) for x in B],
        "C_prefix_0_7": [list(x) for x in C],
        "power_sums_0_7": [list(x) for x in s],
        "scope": "Formal local prefix; no claim that P has all roots in the annulus or that G,Q extend to a global counterexample.",
    }
    Path(__file__).with_name("s3_local_prefix_output.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
