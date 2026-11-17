#!/usr/bin/env python3
"""Independent exact verifier for one p=2, Q=Phi_6 local-lemma counterexample.

This file deliberately does not import the search program.
"""
from __future__ import annotations

R = [1, 0, -1, 0, 1]  # Phi_6(x^2) = 1 - x^2 + x^4
E = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
D = [1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, 0, 0, 1, 1, 1, 1]

def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)

def support(a):
    return [i for i, x in enumerate(a) if x]

def main():
    u = mul(R, E)
    v = mul(R, D)
    h = mul(R, mul(E, D))
    uv = mul(u, v)
    assert u == [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]
    assert v == [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
    assert all(x in (0, 1) for x in u + v + h + uv)
    assert uv == mul(R, h)
    assert E[13] == D[13] == 1 and 13 % 2 == 1
    sums = [i+j for i in support(u) for j in support(v)]
    assert len(sums) == len(set(sums)) == sum(uv)
    print("R support", support(R))
    print("U support", support(u))
    print("V support", support(v))
    print("H support", support(h))
    print("UV support", support(uv))
    print("off-p overlap n=13: E_13=D_13=1")

if __name__ == "__main__":
    main()
