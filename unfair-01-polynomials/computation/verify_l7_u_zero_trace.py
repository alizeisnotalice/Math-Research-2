#!/usr/bin/env python3
"""Exact verification of the first L=7 U-branch zero-propagation fork."""
import json
from pathlib import Path

from sympy import Poly, expand, symbols

HERE = Path(__file__).resolve().parent
trace = json.loads((HERE / "results/021-l7_u_zero_trace.json").read_text())
x, t, u, v = symbols("x t u v")
J = 1 + x + x**2 + x**3 + x**6
C = expand(J * J.subs(x, x**7))
support = [n for n in range(49) if Poly(C, x).coeff_monomial(x**n)]
assert support == trace["target_support"]
assert Poly(C, x).degree() == trace["target_degree"]

a = symbols("a0:20")
b = symbols("b0:20")
A = sum(a[i]*x**i for i in range(20))
B = sum(b[i]*x**i for i in range(20))

def c(n):
    return expand(Poly(expand(A*B), x).coeff_monomial(x**n))

subs = {
    a[0]:1,a[1]:1,a[2]:1,a[3]:1,a[4]:0,a[5]:0,a[6]:t,a[7]:u,a[8]:v,
    b[0]:1,b[1]:0,b[2]:0,b[3]:0,b[4]:0,b[5]:0,b[6]:1-t,b[7]:t-u,b[8]:u-v,
    a[11]:0,b[11]:0,a[12]:0,b[12]:0,
}
subs |= {a[n]: 0 for n in range(20) if n not in support}
subs |= {b[n]: 0 for n in range(20) if n not in support}
assert c(11).subs(subs) == u - v + b[9] + b[10]
step11 = subs | {v:u,b[9]:0,b[10]:0}
assert c(12).subs(step11) == t*(1-t)
step_t0 = step11 | {t:0,u:0,v:0,a[9]:0,b[11]:0,b[12]:0}
assert c(10).subs(step_t0) == a[10]
step_t0 |= {a[10]: 1}
assert c(13).subs(step_t0) == a[13] + b[13]
assert c(19).subs(step_t0) == a[13] + b[16] + b[17]
step_t0 |= {a[13]: 0, b[13]: 1}
assert c(14).subs(step_t0) == a[14] + b[14] + 1
step_t0 |= {a[14]: 0, b[14]: 0}
assert c(15).subs(step_t0) == a[15] + b[15] + 1
step_t0 |= {a[15]: 0, b[15]: 0}
assert c(16).subs(step_t0) == a[16] + b[16] + 2
assert Poly(C, x).coeff_monomial(x**16) == 1
print(json.dumps({"status":"verified","c11":"u-v+b9+b10","c12":"t*(1-t)","t0_c10":"a10","t0_c13":"a13+b13","t0_c16":"a16+b16+2","target_c16":1}, ensure_ascii=False))
