#!/usr/bin/env python3
"""Exact checks for the a=t, b=t^3 sharpness family in (T9a)."""
from fractions import Fraction


def mu(a: int, b: int) -> Fraction:
    return Fraction(a + b - 1, b)


def second(a: int, b: int) -> Fraction:
    return Fraction(1) + Fraction(3 * (a - 1), b) + Fraction(
        2 * (a - 1) * (a - 2), b * (b + 1)
    )


def fibre(a: int, b: int) -> tuple[Fraction, Fraction]:
    m = (mu(a, b) + mu(b, a)) / 2
    s = (second(a, b) + second(b, a)) / 2
    return m, s


# T7 holds on a broad finite exact grid; this is a check, not its proof.
for a in range(1, 101):
    for b in range(1, 101):
        m, s = fibre(a, b)
        assert s <= 4 * m * m

rows = []
for t in (2, 10, 100, 1000):
    m, s = fibre(t, t**3)
    ratio = s / (m * m)
    rows.append((t, str(m), str(s), str(ratio), float(ratio)))
    assert ratio < 4
    if t >= 10:
        assert ratio > Fraction(3, 1)

print("exact grid T7 verification: 1 <= a,b <= 100")
for t, m, s, ratio, approx in rows:
    print(f"t={t}: m={m}, s={s}, s/m^2={ratio} ~= {approx:.12f}")
