#!/usr/bin/env python3
"""Exact Fraction verification of the N=5 r=1,t=5 local-prefix witness."""
from fractions import Fraction as F


def conv(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def main():
    a = [F(1), F(9, 10), F(0), F(181, 1000), F(2421, 40000), F(0)]
    b = [F(1), F(1, 10), F(91, 100), F(0), F(7371, 8000), F(0)]
    c = conv(a, b)
    assert c[:6] == [F(1)] * 6
    assert c[6] == F(220311, 4000000)
    assert 0 < c[6] < 1
    assert a[1] > 0 and b[1] > 0 and a[5] == b[5] == 0
    print("PASS: C_0..C_5 are all 1; C_6=220311/4000000 lies strictly in (0,1)")
    print("scope: exact local-prefix witness only; it is not an original polynomial-factor counterexample")


if __name__ == "__main__":
    main()
