#!/usr/bin/env python3
"""Exact symbolic verification of the L=7 q6 carry-elimination identity."""
import json
from pathlib import Path

from sympy import Poly, expand, symbols

HERE = Path(__file__).resolve().parent
trace = json.loads((HERE / "results/021-l7_carry_elimination_trace.json").read_text())
x, lam, q0, q6 = symbols("x lambda q0 q6")
J = 1 + x + x**2 + x**3 + x**6
P = (lam-q0)*J - q6*x**6
Q = q0 + q6*x**6
remainder = expand(P + J*Q - lam*J)
assert expand(remainder - q6*(x**7 + x**8 + x**9 + x**12)) == 0
C = expand(J*J.subs(x, x**7))
for k in range(1, 8):
    assert Poly(C, x).coeff_monomial(x**((k+1)*7+5)) == 0
print(json.dumps({
    "status": "verified",
    "remainder": str(remainder),
    "zero_target_degrees_checked": [(k+1)*7+5 for k in range(1,8)],
    "result": trace["result"],
}, ensure_ascii=False, indent=2))
