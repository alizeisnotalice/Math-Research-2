#!/usr/bin/env python3
"""Independent exact checks for the Phi_5-only nonuniform allocation note."""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import sympy as sp


def convolve(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, z in enumerate(right):
            out[i + j] += x * z
    return out


def main() -> None:
    y, a = sp.symbols('y a')
    phi5 = [1, 1, 1, 1, 1]
    rows = []
    for k in range(1, 6):
        support_expected = {0}
        for j in range(k):
            support_expected = {old + d * 6**j for old in support_expected for d in range(5)}
        coeff = [int(i in support_expected) for i in range(max(support_expected) + 1)]
        assert all(c in (0, 1) for c in coeff)
        assert {i for i, c in enumerate(coeff) if c} == support_expected
        K = sp.Poly(sum(c * y**i for i, c in enumerate(coeff)), y, domain=sp.ZZ)
        quotient, remainder = sp.div(K, sp.Poly(sum(y**i for i in range(5))**k, y, domain=sp.ZZ))
        assert remainder == 0
        rows.append({'k': k, 'degree': K.degree(), 'support': len(support_expected), 'quotient_degree': quotient.degree()})

    # Check the first-row algebra directly in Q(a), a^2+a-1=0.
    rel = sp.Poly(a*a + a - 1, a)
    for k in range(1, 51):
        alpha = (sp.sqrt(5) - 1) / 2
        m = int(sp.floor(k * alpha))
        # u=v=0 template: E=(1+y)^m and D=Phi_6^(k+m+1).
        e, d = m, -(k + m + 1)
        fp = -e + k*a
        gp = -d - k - k*a
        fp_red = sp.rem(sp.Poly(fp, a), rel).as_expr()
        gp_red = sp.rem(sp.Poly(gp, a), rel).as_expr()
        frac = sp.simplify(k*alpha - m)
        assert sp.simplify(fp_red.subs(a, alpha) - frac) == 0
        assert sp.simplify(gp_red.subs(a, alpha) - (1 - frac)) == 0
        assert e + d + k == -1

    pm = 1-a*y+y*y
    pp = 1+(1+a)*y+y*y
    mm = 1-(1+a)*y+y*y
    mp = 1+a*y+y*y
    def red(expr: sp.Expr) -> sp.Expr:
        return sp.rem(sp.Poly(sp.expand(expr), a), rel).as_expr()
    assert red(pm*pp - sum(y**i for i in range(5))) == 0
    assert red(mm*mp - (y**4-y**3+y**2-y+1)) == 0
    # The Phi5*Phi10 cancellation allocates Phi10 nonuniformly as well.
    R, T = red(pp*mm), red(pm*mp)
    assert sp.Poly(R, y).nth(1) == 0 and sp.Poly(T, y).nth(1) == 0
    assert red(R*T - (1+y**2+y**4+y**6+y**8)) == 0

    here = Path(__file__).resolve()
    report = {
        'status': 'passed_independent_phi5_general_k_replay',
        'K_rows': rows,
        'first_row_k_range': [1, 50],
        'phi5_phi10_cancellation': 'verified; it splits Phi_10 blocks too, hence is outside the Phi_5-only hypothesis',
        'scope': 'Exact integer and Q(sqrt(5)) checks.  The all-k proof is supplied by the base-6 uniqueness and multiplicity argument in phi5_audit.md.',
        'python': sys.version,
        'sympy': sp.__version__,
        'platform': platform.platform(),
        'source_sha256': hashlib.sha256(here.read_bytes()).hexdigest(),
    }
    here.with_suffix('.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'K_rows': rows, 'first_row_k_range': report['first_row_k_range']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
