#!/Users/zzh/.venvs/main/bin/python
"""SymPy audit of the hard-coded H=P*Pbar coefficient formula."""
from pathlib import Path
import hashlib
import sympy as sp

x, a, b, T, N, ell, d, e = sp.symbols("x a b T N ell d e")
abar = T - a
P = 1 + b*x + a*x**2 + (b*a + ell)*x**3 + d*x**4 + e*x**5 + x**6
Pbar = 1 + b*x + abar*x**2 + (b*abar + ell)*x**3 + d*x**4 + e*x**5 + x**6
relation = sp.Poly(a**2 - T*a + N, a)
product = sp.Poly(sp.expand(P*Pbar), x)
reduced = []
for k in range(13):
    rem = sp.Poly(product.coeff_monomial(x**k), a).rem(relation).as_expr()
    assert not rem.has(a), (k, rem)
    reduced.append(sp.expand(rem))

expected = [
    1,
    2*b,
    b*b + T,
    2*b*T + 2*ell,
    2*d + b*b*T + 2*b*ell + N,
    2*e + 2*b*d + 2*b*N + ell*T,
    2 + 2*b*e + d*T + b*b*N + b*ell*T + ell*ell,
    2*b + e*T + b*d*T + 2*d*ell,
    T + e*b*T + 2*e*ell + d*d,
    b*T + 2*ell + 2*d*e,
    2*d + e*e,
    2*e,
    1,
]
assert all(sp.expand(lhs-rhs) == 0 for lhs, rhs in zip(reduced, expected))
out = Path(__file__).with_suffix(".out")
out.write_text(
    "status=PASS\n"
    + "H_ascending=" + repr([str(z) for z in reduced]) + "\n"
    + "script_sha256=" + hashlib.sha256(Path(__file__).read_bytes()).hexdigest() + "\n"
)
print(out.read_text(), end="")
