#!/usr/bin/env python3
"""Read-only, pure-Fraction audit of the parent N=8 certificate."""
from fractions import Fraction as F
from pathlib import Path
import json

CERT = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round044-parent-counterexample/certificate.json")


def add(a, b):
    return [(a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0))
            for i in range(max(len(a), len(b)))]


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def rem(a, f):
    a = list(a)
    while len(a) >= len(f):
        q, shift = a[-1] / f[-1], len(a) - len(f)
        for j, x in enumerate(f):
            a[shift + j] -= q * x
        a.pop()
    return a


def ev(a, x):
    return sum(c * x**i for i, c in enumerate(a))


def main():
    data = json.loads(CERT.read_text())
    assert data["status"] == "passed_exact_local_N8_counterexample"
    assert data["t_minimal_polynomial_low_to_high"] == [216889384781601, -614524809904128, 281474976710656]
    f = [F(x) for x in data["t_minimal_polynomial_low_to_high"]]
    lo, hi = map(F, data["t_isolating_interval"])
    assert lo < hi and ev(f, lo) > 0 > ev(f, hi)
    assert f[1] + 2 * f[2] * hi < 0  # strictly decreasing throughout [lo,hi]

    A = [[F(x) for x in p] for p in data["A_polynomials_in_t"]]
    B = [[F(x) for x in p] for p in data["B_polynomials_in_t"]]
    assert len(A) == len(B) == 9
    assert A[0] == B[0] == [F(1)] and A[8] == B[8] == [F(0)]
    assert A[1] == [F(5, 64)] and B[1] == [F(59, 64)]

    # Recompute every displayed interval, and use it to prove coefficient positivity.
    checked_intervals = []
    for k, p in enumerate(A[:8]):
        assert len(p) <= 2
        left, right = ev(p, lo), ev(p, hi)
        L, U = min(left, right), max(left, right)
        stored = tuple(map(F, data["A_intervals"][k]))
        assert (L, U) == stored
        if k == 0:
            assert L == U == 1
        else:
            assert 0 < L <= U < 1
        checked_intervals.append((L, U))
    assert lo > 0 and hi < 1

    C = [[F(0)] for _ in range(17)]
    for i, ai in enumerate(A):
        for j, bj in enumerate(B):
            C[i + j] = add(C[i + j], mul(ai, bj))
    for k in range(9):
        assert all(x == 0 for x in rem(add(C[k], [-F(1)]), f)), k
    assert C[9] == [F(0)] + A[5]
    # Since t and A5 lie in their certified positive subunit intervals, C9=t*A5 is in (0,1).
    assert 0 < lo * checked_intervals[5][0] < hi * checked_intervals[5][1] < 1
    print("PASS: exact Fraction audit: root isolation, positivity, C[0..8]=1, C9 in (0,1)")


if __name__ == "__main__":
    main()
