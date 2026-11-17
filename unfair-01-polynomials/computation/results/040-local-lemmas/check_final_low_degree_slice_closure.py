#!/usr/bin/env python3
"""Exact regression checks for final_low_degree_slice_closure.md."""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import sys
from pathlib import Path

import sympy as sp


def main() -> None:
    # Complete-cross binary UV permits exactly these low endpoint choices.
    up_vp = [(u, v) for u, v in itertools.product((0, 1), repeat=2) if u + v <= 1]
    assert up_vp == [(0, 0), (0, 1), (1, 0)]

    # Final cubic multiplication, with X standing for x^p.
    X = sp.symbols("X")
    h0, h1, h2, h3 = sp.symbols("h0 h1 h2 h3", nonnegative=True)
    uv = sp.expand((1 + X**2 + X**3) * (h0 + h1 * X + h2 * X**2 + h3 * X**3))
    assert uv.coeff(X, 1) == h1
    assert uv.coeff(X, 3) == h3 + h1 + h0
    assert sp.expand(uv.coeff(X, 3).subs({h0: 1, h1: 1}) - (h3 + 2)) == 0

    # Exact remaining cubic and its unique positive root interval.
    y, a = sp.symbols("y a")
    qcubic = 1 + y**2 + y**3
    assert sp.factor(qcubic) == qcubic
    assert sp.Poly(a**3 + a - 1, a).intervals() == [((sp.Rational(0), sp.Rational(1)), 1)]

    # Phi_5's two real quadratic factors, including the forbidden coefficient.
    phi5 = sum(y**i for i in range(5))
    left = 1 + (1 - sp.sqrt(5)) * y / 2 + y**2
    right = 1 + (1 + sp.sqrt(5)) * y / 2 + y**2
    assert sp.expand(left * right - phi5) == 0
    assert sp.simplify((1 + sp.sqrt(5)) / 2 - 1) > 0

    here = Path(__file__).resolve()
    report = {
        "status": "passed_exact_final_low_degree_closure_regressions",
        "permitted_endpoint_pairs": up_vp,
        "terminal_cubic_identity": "[X^3]((1+X^2+X^3)H)=H_3+H_1+H_0",
        "surviving_before_terminal_step": "Q=1+y^2+y^3, a^3+a-1=0 with a in (0,1)",
        "phi5_factorization_checked": True,
        "scope": "Exact symbolic regressions for the terminal cubic and Phi_5 steps; the universal proof and hypotheses are in the companion manuscript.",
        "python": sys.version,
        "sympy": sp.__version__,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
    }
    here.with_suffix(".json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
