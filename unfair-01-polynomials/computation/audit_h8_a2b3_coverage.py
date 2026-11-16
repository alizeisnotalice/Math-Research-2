#!/usr/bin/env python3
"""Exact combinatorial audit for the H8(x^2)H8(x^3) root-layout cover.

This checks only integer support, local balance solutions, complement orbits,
and the resulting 1482 count.  Galois hypotheses remain analytic inputs.
"""
from __future__ import annotations

from collections import Counter
from itertools import product
import json
from math import comb
from pathlib import Path

H = (1, 1, 1, 0, 0, 0, 0, 1, 1)


def stretch(poly: tuple[int, ...], m: int) -> tuple[int, ...]:
    out = [0] * ((len(poly) - 1) * m + 1)
    out[::m] = poly
    return tuple(out)


def convolve(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return tuple(out)


def main() -> None:
    c = convolve(stretch(H, 2), stretch(H, 3))
    assert set(c) <= {0, 1}
    assert H != tuple(reversed(H))

    # A block has u selected x^2-lift conjugate pairs and v selected
    # x^3-lift conjugate pairs.  Its root-product exponent is E=3u+2v.
    local = [(u, v, 3 * u + 2 * v) for u in range(3) for v in range(4)]
    by_e: dict[int, list[tuple[int, int]]] = {}
    for u, v, e in local:
        by_e.setdefault(e, []).append((u, v))

    assert by_e == {
        0: [(0, 0)], 2: [(0, 1)], 3: [(1, 0)], 4: [(0, 2)],
        5: [(1, 1)], 6: [(0, 3), (2, 0)], 7: [(1, 2)],
        8: [(2, 1)], 9: [(1, 3)], 10: [(2, 2)], 12: [(2, 3)],
    }

    # Complement in the 2+3 groups sends (u,v) to (2-u,3-v), E to 12-E.
    comp = lambda uv: (2 - uv[0], 3 - uv[1])
    assert all(3 * comp(uv)[0] + 2 * comp(uv)[1] == 12 - e
               for e, choices in by_e.items() for uv in choices)

    mixed_by_e = {e: choices for e, choices in by_e.items() if 0 < e < 6}
    # A pair (u,v) labels its multiplicity binom(2,u)binom(3,v), not merely
    # one local choice.  For E<6 the Diophantine solution is unique.
    per_e_counts = {
        e: sum(comb(2, u) * comb(3, v) for u, v in choices) ** 4
        for e, choices in mixed_by_e.items()
    }
    assert per_e_counts == {2: 81, 3: 16, 4: 81, 5: 1296}
    mixed_total = sum(per_e_counts.values())

    terminal = list(product(by_e[6], repeat=4))
    terminal_orbits = {
        min(layout, tuple(comp(uv) for uv in layout)) for layout in terminal
    }
    assert len(terminal) == 16 and len(terminal_orbits) == 8

    f4_2 = sum(comb(2, r) ** 4 for r in range(3))
    f4_3 = sum(comb(3, r) ** 4 for r in range(4))
    total = mixed_total + len(terminal_orbits)
    assert (f4_2, f4_3, total) == (18, 164, 1482)
    assert total == f4_2 * f4_3 // 2 + 2 ** 3 - 2

    report = {
        "status": "passed_exact_combinatorial_coverage_audit",
        "H8_nonreciprocal": True,
        "C_degree": len(c) - 1,
        "C_binary": True,
        "local_balance_solutions": {str(e): choices for e, choices in by_e.items()},
        "mixed_counts_E_2_to_5": per_e_counts,
        "terminal_E6_ordered": len(terminal),
        "terminal_E6_complement_orbits": len(terminal_orbits),
        "franel_F4_2": f4_2,
        "franel_F4_3": f4_3,
        "canonical_nonconstant_complement_representatives": total,
        "scope": "Exact combinatorial audit only; Galois/balance hypotheses are separately proved inputs.",
    }
    Path("audit_h8_a2b3_coverage.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
