#!/usr/bin/env python3
"""Exact finite diagnostics for the proposed factorial-square prefix bound.

This script reads the already certified n=4,5 integer prefix numerators.
For each d it computes the least real C_d satisfying
    p_d <= C_d**d/(d!)**2.
It also records the elementary corner-path lower bound p_{2n-1} >=
4*binom(2n-2,n-1)/(n^2)!; this is deliberately only a non-obstruction
check, not an estimate of p_{2n-1}.
"""
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "reference"


def candidate_constant(p: Fraction, d: int) -> float:
    return float(p.numerator / p.denominator * math.factorial(d) ** 2) ** (1.0 / d)


def one_n(n: int):
    raw = json.loads((SOURCE / f"geodesic-cpp-n{n}.json").read_text())
    a1 = int(raw["A1"])
    assert raw["n"] == n and a1 > 0
    rows = []
    for row in raw["rows"]:
        d, numerator = int(row["d"]), int(row["numerator"])
        p = Fraction(numerator, a1)
        assert 0 < p <= 1
        rows.append({
            "d": d,
            "p": str(p),
            "p_decimal": float(p),
            "required_C": candidate_constant(p, d),
        })
    assert rows[1]["p"] == "1" and rows[2]["p"] == "1/2"
    return {"n": n, "A1": str(a1), "rows": rows,
            "finite_maximum_required_C": max(x["required_C"] for x in rows)}


def corner_lower(n: int):
    # Four orientations from a corner to its opposite corner, and every one
    # has at least one legal completion.  A1<=N! turns this into a lower bound.
    N, d = n*n, 2*n-1
    lower = Fraction(4 * math.comb(2*n-2, n-1), math.factorial(N))
    # RHS with C=e is a convenient scale comparison, not a proposed theorem.
    rhs_e = math.e ** d / math.factorial(d) ** 2
    return {"n": n, "d": d, "p_lower": str(lower),
            "log_lower": math.log(lower.numerator) - math.log(lower.denominator),
            "log_candidate_rhs_C_equals_e": math.log(rhs_e)}


def main():
    out = {"purpose": "finite exact diagnostic; not a uniform bound",
           "definition": "p_d is the ordered geometric-prefix probability under the uniform one-peak law",
           "finite": [one_n(4), one_n(5)],
           "corner_lower_checks": [corner_lower(n) for n in range(4, 31)]}
    (HERE / "candidate-diagnostic.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
