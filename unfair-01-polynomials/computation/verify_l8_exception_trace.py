#!/usr/bin/env python3
"""Independent exact verification of l8_exception_trace.json.

The script uses integer polynomial arithmetic only.  It verifies the target
support, every displayed degree-14 convolution pair, and the forced c16 lower
bound; it does not solve a floating-point or bounded-degree relaxation.
"""
import json
from pathlib import Path

from sympy import Poly, expand, symbols

HERE = Path(__file__).resolve().parent
trace = json.loads((HERE / "results/021-l8_exception_trace.json").read_text())
x = symbols("x")
J = 1 + x + x**2 + x**3 + x**6
C = expand(J * J.subs(x, x**8))
support = [n for n in range(55) if Poly(C, x).coeff_monomial(x**n)]
assert support == trace["target_support"], (support, trace["target_support"])
assert Poly(C, x).degree() == trace["target_degree"]

a = symbols("a0:17")
b = symbols("b0:17")
A = sum(a[i] * x**i for i in range(17))
B = sum(b[i] * x**i for i in range(17))

def coefficient(n):
    return expand(Poly(expand(A * B), x).coeff_monomial(x**n))

prefix = {a[int(k)]: v for k, v in trace["premise_prefix"]["A_0_to_7"].items()}
prefix |= {b[int(k)]: v for k, v in trace["premise_prefix"]["B_0_to_7"].items()}
support_zeros = {a[n]: 0 for n in range(17) if n not in support}
support_zeros |= {b[n]: 0 for n in range(17) if n not in support}

# c8=1, c15=0, c12=0 give these individual nonnegative summands as zero.
forced_zero = {a[8]: 0, b[8]: 0, a[9]: 0, a[14]: 0, b[10]: 0}
subs = prefix | support_zeros | forced_zero

assert coefficient(8).subs(prefix) == a[8] + b[8] + 1
assert coefficient(9).subs(subs) == b[9]
assert coefficient(10).subs(subs) == a[10]

all_pairs_14 = [f"a_{i}*b_{14-i}" for i in range(15)]
trace_pairs_14 = trace["steps"][4]["all_pairs"]
assert all_pairs_14 == trace_pairs_14
assert coefficient(14).subs(subs) == b[14]

subs |= {a[10]: 1, b[14]: 1}
assert coefficient(16).subs(subs) == a[16] + b[16] + 2
assert Poly(C, x).coeff_monomial(x**16) == 1
print(json.dumps({
    "status": "verified",
    "target_degree": Poly(C, x).degree(),
    "target_support": support,
    "c8": "a8+b8+1",
    "c10": "a10",
    "c14": "b14",
    "c16": "a16+b16+2",
    "target_c16": 1,
}, ensure_ascii=False, indent=2))
