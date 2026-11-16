#!/usr/bin/env python3
"""Exact checks for the L=7 EA and V non-block branch eliminations."""
import json
from pathlib import Path

from sympy import Poly, Rational, expand, symbols

HERE = Path(__file__).resolve().parent
trace = json.loads((HERE / "results/021-l7_nonblock_elimination_trace.json").read_text())
x, t, u, v = symbols("x t u v")
J = 1 + x + x**2 + x**3 + x**6
C = expand(J * J.subs(x, x**7))
S = [n for n in range(49) if Poly(C, x).coeff_monomial(x**n)]
assert S == trace["support"] and Poly(C, x).degree() == trace["target_degree"]

def setup(branch):
    a = symbols("a0:20")
    b = symbols("b0:20")
    A = sum(a[i]*x**i for i in range(20)); B = sum(b[i]*x**i for i in range(20))
    def c(n): return expand(Poly(expand(A*B), x).coeff_monomial(x**n))
    z = {a[n]:0 for n in range(20) if n not in S} | {b[n]:0 for n in range(20) if n not in S}
    if branch == "V":
        z |= {a[0]:1,a[1]:0,a[2]:1,a[3]:0,a[4]:0,a[5]:0,a[6]:t,a[7]:u,a[8]:v,
              b[0]:1,b[1]:1,b[2]:0,b[3]:0,b[4]:0,b[5]:0,b[6]:1-t,b[7]:1-t-u,b[8]:t-u-v}
    else:
        z |= {a[0]:1,a[1]:1,a[2]:0,a[3]:0,a[4]:0,a[5]:0,a[6]:t,a[7]:u,a[8]:v,
              b[0]:1,b[1]:0,b[2]:1,b[3]:0,b[4]:0,b[5]:0,b[6]:1-t,b[7]:t-u,b[8]:1-2*t+u-v}
    return a,b,c,z

# V: force t endpoints; each endpoint contradicts an exact coefficient identity.
a,b,c,z = setup("V")
assert c(11).subs(z) == a[10] + b[9]
assert c(12).subs(z) == b[10] + t*(1-t)
z0 = z | {t:0,u:0,v:0,a[9]:0,a[10]:0,b[9]:0,b[10]:0}
assert c(9).subs(z | {t:0,u:0,v:0,a[10]:0,b[9]:0,b[10]:0}) == a[9] + 1
assert c(10).subs(z0) == 0
v1 = z | {t:1,u:0,a[10]:0,b[9]:0,b[10]:0}
assert c(9).subs(v1) == a[9] + v
assert c(10).subs(v1) == a[9] - v + 1
z1 = v1 | {v:Rational(1,2),a[9]:Rational(1,2),b[16]:0,b[17]:0,b[13]:0}
assert c(13).subs(z1) == a[13]
assert c(14).subs(z1 | {a[13]:1}) == a[14] + b[14] + Rational(3,2)

# EA: same endpoint split, with the other orientation of half-coefficients.
a,b,c,z = setup("EA")
assert c(11).subs(z) == a[9] + b[10]
assert c(12).subs(z) == a[10] + t*(1-t)
z1 = z | {t:1,u:1,v:0,a[9]:0,a[10]:0,b[9]:0,b[10]:0}
assert c(10).subs(z1) == 0
e0 = z | {t:0,u:0,a[9]:0,a[10]:0,b[10]:0}
assert c(9).subs(e0) == b[9] - v + 1
assert c(10).subs(e0) == b[9] + v
z0 = e0 | {v:Rational(1,2),b[9]:Rational(1,2),a[13]:0,a[17]:0}
assert c(13).subs(z0) == b[13]
assert c(14).subs(z0 | {b[13]:1}) == a[14] + b[14] + Rational(3,2)
print(json.dumps({"status":"verified","branches":["V","EA"],"conclusion":trace["result"]}, ensure_ascii=False))
