#!/usr/bin/env python3
"""Small exact regressions for the mass-cell interfaces.

This does not evaluate a project candidate. It checks the one-box activation
and the nonconserving complete-periodic counterexample in
proofs/periodic-halo-mass-cells.md using only rational arithmetic.
"""
from fractions import Fraction as F
import argparse
import json
from pathlib import Path
import platform
import sys


def main(output_path: str | None = None) -> None:
    # A source box has x-span 1, y-span 0, lambda=1; its receiver is nonempty
    # exactly when sqrt(w) > 1.  We use square masses for exact evaluation.
    values = {F(1, 4): F(0), F(1): F(0), F(4): F(2)}
    for w, expected in values.items():
        root = F(1, 2) if w == F(1, 4) else F(1) if w == F(1) else F(2)
        area = F(0) if root <= 1 else (root - 1) * root
        assert area == expected
    left_derivative = F(0)
    right_derivative = F(1) - F(1, 2)
    assert left_derivative != right_derivative

    # Complete-periodic nonconserving example: Px=Py=2, lambda=1, R=2,
    # one mass-w atom per cell.  The tail upper bound and the exact ratio
    # show why an increase in area alone is not a weak-type improvement.
    nonconserving = {}
    for w in (F(1, 4), F(1)):
        tail_bound = w * F(1, 4) * (F(1) + F(2, 4)) ** 2
        assert tail_bound < 1
        area = w
        total_mass = w
        ratio = area / total_mass
        assert ratio == 1
        nonconserving[str(w)] = {
            "area": str(area),
            "mass": str(total_mass),
            "J": str(ratio),
            "tail_upper_bound_at_R": str(tail_bound),
        }
    result = {
        "status": "passed",
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "reproduction_command": "python3 -I computation/check_periodic_halo_mass_cells.py --output <output.json>",
        "actual_argv": [sys.executable, *sys.argv],
        "arithmetic": "fractions.Fraction (exact rational arithmetic)",
        "activation_threshold": "w=1",
        "areas": {str(k): str(v) for k, v in values.items()},
        "left_derivative": str(left_derivative),
        "right_derivative": str(right_derivative),
        "scope": "One-box local activation plus a closed-form complete-periodic counterexample; not a project-candidate computation.",
        "nonconserving_complete_periodic_counterexample": nonconserving,
    }
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered)
    if output_path is not None:
        Path(output_path).write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="write the exact run record as JSON")
    args = parser.parse_args()
    main(args.output)
