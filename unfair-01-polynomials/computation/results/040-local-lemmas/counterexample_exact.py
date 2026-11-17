#!/Users/zzh/.venvs/main/bin/python
"""Exact certificate for the round-40 direct-sum gcd counterexample."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


x = sp.symbols("x")


def main() -> None:
    A = (0, 3, 7, 11, 14)
    B = (0, 1, 13, 23, 33, 45, 46)
    U = sp.Poly(sum(x**i for i in A), x, domain=sp.ZZ)
    V = sp.Poly(sum(x**i for i in B), x, domain=sp.ZZ)
    sums = sorted(a + b for a in A for b in B)
    assert len(sums) == len(set(sums)) == 35
    UV = U * V
    assert all(c in (0, 1) for c in UV.all_coeffs())

    ell = sp.Poly(x**10 - x**9 + x**7 - x**6 + x**5 - x**4 + x**3 - x + 1, x, domain=sp.ZZ)
    phi = {n: sp.Poly(sp.cyclotomic_poly(n, x), x, domain=sp.ZZ) for n in (5, 6, 7, 30, 66)}
    assert sp.Poly(U.as_expr() - ell.as_expr() * phi[5].as_expr(), x, domain=sp.ZZ).is_zero
    assert sp.Poly(V.as_expr() - ell.as_expr() * sp.prod(phi[n].as_expr() for n in (6, 7, 30, 66)), x, domain=sp.ZZ).is_zero
    assert sp.gcd(U, V) == ell
    assert ell.eval(0) == 1 and ell.eval(-1) == -1

    result = {
        "arithmetic": "exact SymPy ZZ",
        "A": list(A), "B": list(B), "direct_sum_size": len(sums), "UV_support": sums,
        "gcd": str(ell.as_expr()),
        "factorization_U": "ell*Phi_5",
        "factorization_V": "ell*Phi_6*Phi_7*Phi_30*Phi_66",
        "noncyclotomic_witness": "ell(0)=1 and ell(-1)=-1, so ell has a root in (-1,0)",
    }
    path = Path(__file__).with_name("counterexample_exact.json")
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
