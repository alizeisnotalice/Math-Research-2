#!/usr/bin/env python3
"""Exact small-grid audit of p_n=P_1(the top row is the first n-set).

For n<=4 this computes both A_1 by the connected-order recurrence and
H_n by the reverse completion recurrence, with Python integers only.
"""
from fractions import Fraction
import json
from pathlib import Path


def calc(n):
    N = n * n
    full = (1 << N) - 1
    adj = [0] * N
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                if 0 <= X < n and 0 <= Y < n:
                    adj[v] |= 1 << (Y * n + X)

    # a[S] counts connected growth words on the induced set S.
    a = [0] * (full + 1)
    for v in range(N):
        a[1 << v] = 1
    for size in range(2, N + 1):
        for S in range(1, full + 1):
            if S.bit_count() != size:
                continue
            z = S
            total = 0
            while z:
                bit = z & -z
                v = bit.bit_length() - 1
                T = S ^ bit
                if adj[v] & T:
                    total += a[T]
                z ^= bit
            a[S] = total

    # h[S] counts completions from a connected revealed set S.
    h = [0] * (full + 1)
    h[full] = 1
    for size in range(N - 1, 0, -1):
        for S in range(1, full + 1):
            if S.bit_count() != size:
                continue
            boundary = 0
            z = S
            while z:
                bit = z & -z
                boundary |= adj[bit.bit_length() - 1]
                z ^= bit
            boundary &= full ^ S
            total = 0
            z = boundary
            while z:
                bit = z & -z
                total += h[S | bit]
                z ^= bit
            h[S] = total

    top = (1 << n) - 1
    A1 = a[full]
    H = h[top]
    two_rows = (1 << (2 * n)) - 1
    HR = h[two_rows]
    p = Fraction((1 << (n - 1)) * H, A1)
    return {"n": n, "A1": A1, "H": H, "p": str(p),
            "p_float": float(p), "factorial_scaled": str(p * _factorial(n)),
            "factorial_scaled_float": float(p * _factorial(n)),
            "H_two_rows": HR, "H_over_H_two_rows": str(Fraction(H, HR)),
            "two_row_lower_bound": _factorial(n) * HR,
            "two_row_lower_bound_verified": H >= _factorial(n) * HR}


def _factorial(n):
    x = 1
    for k in range(2, n + 1):
        x *= k
    return x


out = {"results": [calc(n) for n in range(2, 5)]}
Path(__file__).with_name("toprow_singlepeak_probability_n2_n4.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
