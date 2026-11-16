#!/usr/bin/env python3
"""Exact-rational audit of the all-n numerical inequalities used in r014."""

from fractions import Fraction
import json
from pathlib import Path
import sys


def B(N, delta=4):
    return sum(
        (Fraction(N - 1 - delta * q, N - q - 1) if N - 1 - delta * q > 0 else 0)
        for q in range(1, N - 1)
    )


def main():
    rows = []
    for n in range(4, 101):
        N = n * n
        b = B(N)
        coarse = Fraction(N, 8) - Fraction(3, 2)
        assert b >= coarse
        rows.append({"n": n, "N": N, "certificate": str(coarse), "B_ge_certificate": True})
    ans = {"range": "n=4,...,100", "all_exact_rational_checks_passed": True,
           "first": rows[0], "last": rows[-1], "status": "passed"}
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("ratio-lower-bound-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
