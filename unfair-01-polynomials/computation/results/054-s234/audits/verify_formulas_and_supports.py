"""Exact symbolic audit for the round-054 S234 draft.

This is an audit aid only: it expands the displayed Newton and low-row
identities and independently enumerates the support reduction used in §6.
"""
from itertools import combinations
from sympy import Symbol, Poly, expand, symbols

x = Symbol("x")
L, a, k, d, e, f = symbols("L a k d e f")
g1, g2, g3, g4, g5, g6, g7 = symbols("g1:8")

P = 1 + L*x + a*x**2 + (L*a+k)*x**3 + d*x**4 + e*x**5 + f*x**6 + x**7
G = 1 + sum(g*x**i for i, g in enumerate((g1,g2,g3,g4,g5,g6,g7), 1))
A = expand(P*G)

subs_low = {g1:-L, g3:-L*g2-k}

def coeff(poly, j):
    return expand(poly).coeff(x, j)

u2 = g2-L**2
u4 = g4-L**2*g2-2*L*k
eta = g5+e+k*g2+L*g4
rho = g6+f-2*L*e+L*eta-2*L*k*g2-L**2*g4-k**2
u6 = a*(u4+L*k)+d*(u2+L**2)
R7 = 1+L*rho-L**2*eta+2*L**2*k*g2+(L**3+k)*g4+L*k**2+e*(2*L**2+g2)-2*L*f
u7 = a*(eta-e-k*g2)-d*(L*g2+k)

checks = {
    "A1_zero": coeff(A,1).subs(subs_low),
    "A3_zero": coeff(A,3).subs(subs_low),
    "A2": expand(coeff(A,2).subs(subs_low)-(a+u2)),
    "A4": expand(coeff(A,4).subs(subs_low)-(d+a*u2+u4)),
    "A5": expand(coeff(A,5).subs(subs_low)-(-k*a-L*d+eta)),
    "A6": expand(coeff(A,6).subs(subs_low)-(u6+rho)),
    "A7": expand(coeff(A,7).subs(subs_low)-(u7+(g7+R7))),
}

# Newton sums for P(x)=prod(1-w_i x).
ps = [None] + [Poly(P, x).coeff_monomial(x**j) for j in range(1,8)]
s = [None]
for j in range(1,8):
    s.append(expand(-sum(ps[i]*s[j-i] for i in range(1,j)) - j*ps[j]))
expected = {
    2:L**2-2*a,
    3:-L**3-3*k,
    4:2*a**2+L**4+4*L*k-4*d,
    5:5*k*a+5*L*d-L**5-5*L**2*k-5*e,
    6:-2*a**3-6*L*k*a+6*a*d-6*L**2*d+L**6+6*L**3*k+3*k**2+6*L*e-6*f,
    7:-7*k*a**2+7*L**2*k*a-7*L*a*d+7*a*e-L**7-7*L**4*k+7*L**3*d-7*L**2*e-7*L*k**2+7*L*f+7*k*d-7,
}
for j in range(2,8):
    checks[f"s{j}"] = expand(s[j]-expected[j])

assert all(v == 0 for v in checks.values()), checks

# Independent support enumeration in draft §6.
U = set(range(1,7))
remaining = []
for n in range(1,7):
    for tpl in combinations(sorted(U), n):
        S = set(tpl)
        r, t = min(S), 7-max(S)
        forbidden = (
            r == t or
            (r == 1 and t in {2,3,4}) or
            (t == 1 and r in {2,3,4}) or
            (r,t) in {(2,4),(4,2)}
        )
        if not forbidden:
            remaining.append(tuple(sorted(S)))
expected_remaining = [
    (1,), (2,), (3,), (4,), (5,), (6,),
    (1,2), (2,4), (3,5), (5,6), (2,3,4), (3,4,5),
]
assert remaining == expected_remaining, (remaining, expected_remaining)

print("sympy", __import__("sympy").__version__)
print("all displayed low-row and Newton identities: exact PASS")
print("support reduction remaining", remaining)
