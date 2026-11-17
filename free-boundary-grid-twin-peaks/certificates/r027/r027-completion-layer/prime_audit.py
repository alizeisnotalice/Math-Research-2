#!/usr/bin/env python3
"""Deterministic 64-bit Miller--Rabin audit and CRT-width calculation."""
from __future__ import annotations
from math import factorial
from pathlib import Path
import json

BASES = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)

def prime64(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while not d & 1:
        s += 1; d //= 2
    for a in BASES:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

# First four primes below 2^63, searched by this script rather than asserted.
found = []
x = (1 << 63) - 1
while len(found) < 4:
    if prime64(x):
        found.append(x)
    x -= 2
P = 1
for p in found: P *= p
out = {
    "method": "deterministic 64-bit Miller-Rabin bases from the proven n<2^64 basis set",
    "bases": list(BASES), "primes": found,
    "all_primes_verified": all(prime64(p) for p in found),
    "minimum_prime": min(found), "all_exceed_56": min(found) > 56,
    "product": str(P), "product_bit_length": P.bit_length(),
    "factorial_56": str(factorial(56)), "factorial_56_bit_length": factorial(56).bit_length(),
    "product_exceeds_factorial_56": P > factorial(56),
}
Path("prime-audit.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
