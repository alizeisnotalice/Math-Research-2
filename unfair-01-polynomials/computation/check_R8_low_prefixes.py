#!/usr/bin/env python3
"""Exact checks for the six classified R8 low 0..8 prefixes."""
from sympy import Poly, expand, symbols

x, s = symbols("x s")
T = sum(x**i for i in range(6)) + x**8
prefixes = [
    (1+x+x**2+x**3+x**4+x**5, 1),
    (1+x+x**2, 1+x**3),
    (1+x, 1+x**2+x**4),
]
for A0, B0 in prefixes:
    for A, B in ((A0+s*x**8, B0+(1-s)*x**8),
                 (B0+(1-s)*x**8, A0+s*x**8)):
        for n in range(9):
            assert Poly(expand(A*B-T), x).coeff_monomial(x**n) == 0
print("verified 6 oriented low-prefix families; s+(1-s)=1 at degree 8")
