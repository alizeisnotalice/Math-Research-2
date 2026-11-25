"""Audit the two-prime-row obstruction at gap ceil(N/4)."""

import json
from math import comb, gcd, isqrt, prod, lcm
from fractions import Fraction
from lattice import column
from intervals import harmonic


def prime(n):
    return n > 1 and all(n % p for p in range(2, isqrt(n) + 1))


def moments(N, m):
    U = 0
    V = Fraction()
    for i in range(m + 1):
        c = comb(N, i) ** 2
        U += c
        V += c * (harmonic(i) + harmonic(N) - harmonic(N - i))
    return U, V


def check(N):
    g = (N + 3) // 4
    a, b = [p for p in range(N + 1, N + N // 2 + 1) if prime(p)][:2]
    U, V = moments(N, a - N - 1)
    W, Z = moments(N, b - N - 1)
    delta = U * Z - W * V
    assert delta > 0
    E, _, _ = column(N)
    d = lcm(*range(1, 2 * N + 1))
    assert E[a] == 2 * d * (U * harmonic(N) - V) and E[b] == 2 * d * (
        W * harmonic(N) - Z
    )
    L = prod(p for p in range(2 * N + 1, 2 * (N + g) + 1) if prime(p))
    bad = gcd(L, E[a], E[b])
    assert delta.numerator % bad == 0 and gcd(L, delta.denominator) == 1
    required = L // gcd(L, *E.values())
    assert required % (L // gcd(L, delta.numerator)) == 0
    if N <= 64:
        for j in range(N + g, 2 * N + 1):
            assert all(v % L == 0 for v in column(j)[0].values())
    return dict(
        N=N,
        gap=g,
        rows=[a, b],
        L=L,
        bad_two=bad,
        necessary_divisor=required,
        det_numerator_bits=delta.numerator.bit_length(),
    )


if __name__ == "__main__":
    print(
        json.dumps(
            [check(n) for n in (32, 40, 48, 64, 96, 128, 256, 512, 1024)], indent=2
        )
    )
