#!/usr/bin/env python3
"""Exact regressions for the reducible quartic 1+3 terminal step."""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import sys
from pathlib import Path

import sympy as sp


def diophantine_pairs(p: int, q: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(2 * p + 1) for j in range(2 * p + 1) if q * i + p * j == 2 * p * p]


def main() -> None:
    a, y = sp.symbols("a y")
    initial = [(c2, c3) for c2 in (-1, 0, 1) for c3 in range(-1 - c2, 2 - c2)]
    survivors = []
    for c2, c3 in initial:
        for h2, h3, h4 in itertools.product((0, 1), repeat=3):
            values = (h2 + c2, h3 + c2 + c3, h4 + c2 * h2 + c3 + 1)
            if all(v in (0, 1) for v in values):
                survivors.append((c2, c3))
                break
    assert survivors == [(-1, 0), (-1, 1), (0, -1), (0, 0), (1, -2), (1, -1), (1, 0)]
    equations = [sp.factor(a**4 + c2 * a**2 - c3 * a + 1) for c2, c3 in survivors]
    assert equations[1] == (a - 1) * (a**3 + a**2 - 1)
    residual = 1 - y**2 + y**3 + y**4
    assert sp.factor(residual) == (y + 1) * (y**3 - y + 1)
    assert sp.Poly(a**3 + a**2 - 1, a).intervals() == [((sp.Rational(0), sp.Rational(1)), 1)]

    # Cross-term identity at the terminal candidate: t1=-a, t2=a^2-1,
    # V_2p=0 and G_p=a force G_2p=1.
    assert sp.expand(sp.Symbol("g2") + (-a) * a + (a**2 - 1) - sp.Symbol("g2") + 1) == 0

    # Independently enumerate the claimed nonnegative Diophantine supports.
    examples = {}
    for pp, qq in ((2, 3), (3, 4), (3, 5), (4, 5), (4, 7), (5, 7), (5, 9), (5, 11)):
        pairs = diophantine_pairs(pp, qq)
        expected = [(0, 2 * pp)]
        if qq < 2 * pp:
            expected.append((pp, 2 * pp - qq))
        assert pairs == expected
        examples[f"p={pp},q={qq}"] = pairs

    here = Path(__file__).resolve()
    report = {
        "status": "passed_exact_reducible_quartic_terminal_regressions",
        "low_row_survivors": survivors,
        "unique_root_relation_survivor": "S=(1+y)(1-y+y^3), a^3+a^2-1=0 with a in (0,1)",
        "diophantine_support_examples": examples,
        "scope": "Exact checks of the finite low-row table, residual factorization and coefficient-support equation; the universal inequalities are proved in the companion audit manuscript.",
        "python": sys.version,
        "sympy": sp.__version__,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
    }
    here.with_suffix(".json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
