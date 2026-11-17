#!/Users/zzh/.venvs/main/bin/python
"""Exact algebraic regressions for the conditional quartic factor-type proof."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import sympy as sp


def main() -> None:
    a, y = sp.symbols("a y")

    # 1+3: the nine pairs from c2 and c2+c3 in {-1,0,1}, then the fourth row.
    initial = [(c2, c3) for c2 in (-1, 0, 1) for c3 in range(-1 - c2, 2 - c2)]
    assert len(initial) == 9
    surviving = []
    for c2, c3 in initial:
        for h2, h3, h4 in itertools.product((0, 1), repeat=3):
            rows = (h2 + c2, h3 + c2 + c3, h4 + c2 * h2 + c3 + 1)
            if all(value in (0, 1) for value in rows):
                surviving.append((c2, c3))
                break
    assert surviving == [(-1, 0), (-1, 1), (0, -1), (0, 0), (1, -2), (1, -1), (1, 0)]
    equations = [sp.factor(a**4 + c2 * a**2 - c3 * a + 1) for c2, c3 in surviving]
    assert equations[1] == (a - 1) * (a**3 + a**2 - 1)
    assert sp.factor(y**4 + y**3 - y**2 + 1) == (y + 1) * (y**3 - y + 1)

    # 2+2: exact terminal factorisations arising from equality in the
    # a+a^{-1} >= 2 and rt >= -1 / 0 budgets.
    assert sp.expand((1 + y + y**2) * (1 - y + y**2)) == 1 + y**2 + y**4
    assert sp.expand((1 - y + y**2) * (1 + y**2)) == 1 - y + 2 * y**2 - y**3 + y**4

    result = {
        "arithmetic": "exact integer enumeration and SymPy identities",
        "one_plus_three_initial_pairs": initial,
        "one_plus_three_after_fourth_row": surviving,
        "one_plus_three_root_equations": [str(item) for item in equations],
        "two_plus_two_terminal_factorizations": {
            "c1_0_c2_1": "(1+y+y^2)(1-y+y^2)=1+y^2+y^4",
            "c1_minus1_c2_2": "(1-y+y^2)(1+y^2)=1-y+2y^2-y^3+y^4",
        },
        "scope": "Checks the displayed low-row algebra only; the low-end boxes and integer-slice descent are external hypotheses.",
    }
    path = Path(__file__).with_name("audit_quartic_factor_types.json")
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
