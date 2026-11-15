#!/usr/bin/env python3
"""Exact minimal counterexample to the unqualified shape-free R66 condition."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import time


def main():
    start = time.perf_counter()
    xs = [F(0), F(1, 4), F(1, 2), F(3, 4), F(67, 64)]
    masses = [F(1, 4)] * 4 + [F(1, 64)]
    rows = []
    for mask in range(1, 32):
        ids = [i for i in range(5) if mask >> i & 1]
        mass = sum(masses[i] for i in ids)
        span = max(xs[i] for i in ids) - min(xs[i] for i in ids)
        closure = sum(1 << i for i, x in enumerate(xs)
                      if min(xs[j] for j in ids) <= x <= max(xs[j] for j in ids))
        rows.append({"mask": mask, "mass": str(mass), "span_squared": str(span * span),
                     "eligible": span * span < mass, "closed": closure == mask})
    eligible = [r for r in rows if r["eligible"]]
    L = max(F(r["mass"]) for r in eligible)
    maximizing = [r["mask"] for r in eligible if F(r["mass"]) == L]
    high = [r["mask"] for r in eligible if r["closed"] and F(r["mass"]) > F(3, 4) * L and r["mask"] & 15]
    assert L == 1 and maximizing == [15] and high == [15, 30]
    reached = 1
    while True:
        nxt = reached
        for r in eligible:
            if r["mask"] & nxt:
                nxt |= r["mask"]
        if nxt == reached:
            break
        reached = nxt
    assert reached == 31
    total = sum(masses)
    nonterminal = (xs[-1] - xs[0]) ** 2 - F(35, 33) * total
    assert nonterminal > 0
    assert 3 * 428 ** 2 - 741 ** 2 == 471 > 0
    assert F(3) < F(49, 16)
    assert F(539, 512) + F(27, 128) * F(7, 4) == F(91, 64) < F(3, 2)
    assert F(3) > F(23, 16) ** 2  # remaining halo right edge lies inside C
    data = {
        "schema": "r066-five-source-shape-loss-v1",
        "status": "PASS exact source enumeration and scalar comparisons",
        "scope": "Refutes shape-free compensation condition on all connected nonterminal inputs; S is thin, so not a fat-branch counterexample or a 3.5 deletion counterexample.",
        "points": [[str(x), "0"] for x in xs],
        "initial_masses": [str(q) for q in masses],
        "masses": [str(q) for q in masses],
        "S_indices": [0, 1, 2, 3], "S_mask": 15, "theta": "3/4",
        "all_nonempty_subsets": rows,
        "maximum_eligible_mass": str(L), "maximum_masks": maximizing,
        "high_closed_masks": high, "eligible_connected": True,
        "actual_total_mass": str(total), "actual_bbox_span_squared": str(xs[-1] ** 2),
        "nonterminal_margin": str(nonterminal), "S_widths": ["3/4", "0"],
        "fat_branch": False,
        "quantities_as_rational_plus_coefficient_sqrt3": {
            "C_area": ["3/4", "3/8"],
            "low_old_credit": ["1/64", "0"],
            "high_unpaid": ["163/512", "-21/128"],
            "compensated_high": ["155/512", "-21/128"],
            "shape_free_margin": ["7/4", "-1"],
            "compensated_high_minus_shape_free_margin": ["-741/512", "107/128"],
            "actual_shape_budget": ["539/512", "27/128"],
        },
        "positive_comparison_integer": 471,
        "actual_shape_budget_strict_upper": "91/64",
        "python": sys.version, "arithmetic": "integer/Fraction; no floats or randomness",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "seconds": time.perf_counter() - start,
    }
    out = Path(sys.argv[1]);out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({"subsets": 31, "max": maximizing, "high": high,
                      "nonterminal_margin": str(nonterminal), "positive_integer": 471,
                      "seconds": data["seconds"]}))
    print(hashlib.sha256(out.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
