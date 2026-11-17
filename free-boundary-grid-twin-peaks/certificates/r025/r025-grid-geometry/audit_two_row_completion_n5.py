#!/usr/bin/env python3
"""Exact completion count from the first two rows of P_5 square.

The state is only the 15 unrevealed vertices, so this is a 2^15 integer DP.
"""
import json
from math import factorial
from pathlib import Path

n = 5
N = n * n
initial = (1 << (2 * n)) - 1
full = (1 << N) - 1
adj = [0] * N
for y in range(n):
    for x in range(n):
        v = y * n + x
        for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
            if 0 <= X < n and 0 <= Y < n:
                adj[v] |= 1 << (Y * n + X)

U = full ^ initial
dp = {0: 1}
for size in range(1, U.bit_count() + 1):
    S = U
    while True:
        if S.bit_count() != size:
            if S == 0:
                break
            S = (S - 1) & U
            continue
        total = 0
        z = S
        while z:
            bit = z & -z
            v = bit.bit_length() - 1
            prior = initial | (S ^ bit)
            if adj[v] & prior:
                total += dp[S ^ bit]
            z ^= bit
        dp[S] = total
        if S == 0:
            break
        S = (S - 1) & U

H_R = dp[U]
H_T = 153450311313600  # independently exact r024 top-row completion count
out = {
    "n": n,
    "H_first_two_rows": H_R,
    "H_top_row_imported_exact": H_T,
    "factorial_n": factorial(n),
    "normalized_c_n": {"numerator": H_T, "denominator": factorial(n) * H_R},
    "lower_bound_verified": H_T >= factorial(n) * H_R,
}
Path(__file__).with_name("two_row_completion_n5.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
