"""Exact support and low-row convolution audit for H_8(x^3)H_8(x^5).

No optimization or numerical search is used.  The output lists all nine
forward aliases, derives the four transport rows, and prints the exact
coefficient equations for x_0*y_0 on the allowed row-zero support.
"""
from collections import Counter, defaultdict
import sympy as sp

d = 8

def p(a):
    return Counter({(0, 0): 1, (0, a): 1, (0, 2*a): 1,
                    (a, 0): 1, (a, -a): 1})

def mul(A, B):
    C = Counter()
    for (r, e), c in A.items():
        for (s, f), q in B.items():
            C[r+s, e+f] += c*q
    return C

Tcnt = mul(p(3), p(5))
T = set(Tcnt)
im = defaultdict(list)
for r, e in T:
    im[r*d+e].append((r, e))
assert all(len(v) == 1 for v in im.values())
actual = mul(Tcnt, Tcnt)
aliases = [
    ((r, e), im[r*d+e][0])
    for r, e in sorted(actual)
    if r*d+e in im and im[r*d+e][0] != (r, e)
]

by_pair = defaultdict(list)
for source, target in aliases:
    by_pair[(source[0], target[0])].append((source[1], target[1]))
print("target rows:", {r: sorted(e for rr, e in T if rr == r)
                       for r in sorted({r for r, _ in T})})
print("aliases by layer:", {pair: coords for pair, coords in by_pair.items()})

# Correct expected values directly from the target coordinates.
assert set(by_pair) == {(0, 3), (3, 5), (5, 8)}
assert {pair: len(coords) for pair, coords in by_pair.items()} == {
    (0, 3): 4, (3, 5): 3, (5, 8): 2
}

# Formal four-transport row bookkeeping.
L3, L5, U3, U5, F, D, E, z = sp.symbols("L3 L5 U3 U5 F D E z")
V0 = L3*L5 + z**3*F
R3 = U3*L5 - F + z**2*D
S5 = L3*U5 - D + z**3*E
W8 = U3*U5 - E
print("row template:",
      {"K0": V0, "K3": R3, "K5": S5, "K8": W8, "z": "u^8"})

# In the non-pure branch K6=0 forces x3*y3=0, and W!=0 gives the
# rank-one identity K0*K8=K3*K5.  After the middle transport D vanishes,
# its exact residual factors as follows (F is the 0->3 transport and E is
# the 5->8 transport).
rank_residual_E0 = sp.factor(
    (L3*L5 + z**3*F) * (U3*U5 - E)
    - (U3*L5 - F) * (L3*U5 + z**3*E)
)
assert sp.expand(
    rank_residual_E0 + (L3 + U3*z**3) * (L5*E - U5*F)
) == 0
print("non-pure rank factor after middle transport vanishes:", rank_residual_E0)

# Exact row-zero coefficient formulas.
S0 = (0, 3, 5, 6, 8, 10, 11, 13, 16)
a = {i: sp.symbols(f"a{i}") for i in S0}
b = {i: sp.symbols(f"b{i}") for i in S0}
def coeff(n):
    return sp.expand(sum(a[i]*b[n-i] for i in S0 if n-i in b))
for n in range(0, 33):
    print(f"K0[{n}] = {coeff(n)}")

# First-fraction obstruction at degree 3:
# Let a3=t, b3=1-t.  K0[9]=a3*b6+a6*b3=0.
t = sp.symbols("t", real=True)
eq6 = sp.expand(coeff(6).subs({a[3]: t, b[3]: 1-t}))
eq9 = sp.expand(coeff(9).subs({a[3]: t, b[3]: 1-t}))
print("after a3=t,b3=1-t: K0[6] =", eq6, "; K0[9] =", eq9)
print("If 0<t<1, K0[9]=0 forces a6=b6=0, contradicting K0[6]=1.")

# Pure-v^0 branch after the S_8 root-balance reduction:
# q is monic degree eight and divides u^8 W, hence q=1+b*u^3+u^5+u^8.
# Exact polynomial division of the full target has a nonzero resultant
# already modulo 3 between the constant and x coefficients of its remainder.
x, beta = sp.symbols("x beta")
H8 = lambda z0: 1 + z0 + z0**2 + z0**7 + z0**8
C8 = sp.expand(H8(x**3) * H8(x**5))
qbeta = 1 + beta*x**3 + x**5 + x**8
remainder = sp.Poly(sp.div(C8, qbeta, x)[1], x)
r0, r1 = remainder.nth(0), remainder.nth(1)
assert sp.Poly(sp.resultant(r0, r1, beta), beta).is_zero is False
assert int(sp.resultant(r0, r1, beta)) % 3 == 2
assert sp.gcd(sp.Poly(r0, beta, modulus=3), sp.Poly(r1, beta, modulus=3)).degree() == 0
print("pure degree-8 candidate q_beta remainder resultant mod 3 = 2; no beta can divide C8.")
