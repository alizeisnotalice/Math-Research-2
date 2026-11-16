#!/usr/bin/env python3
"""Relative-path archival replay of the r022 exact geodesic recurrence."""
from fractions import Fraction
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"


def rows_for(n):
    candidate = INPUTS / f"geodesic-cpp-n{n}.json"
    if candidate.exists():
        data = json.loads(candidate.read_text())
        return int(data["A1"]), data["rows"]
    data = json.loads((INPUTS / "geodesic-r001.json").read_text())
    row = next(x for x in data["results"] if x["n"] == n)
    return int(row["A1"]), row["rows"]


def one(n):
    a1, rows = rows_for(n)
    p = {int(r["d"]): Fraction(int(r["numerator"]), a1) for r in rows}
    kernel = json.loads((INPUTS / f"n{n}-extension-analysis.json").read_text())
    q_kernel = {int(r["prefix_vertices"]): Fraction(r["geodesic_extension_probability"])
                for r in kernel["rows"]}
    checks = []
    product = Fraction(1)
    for d in range(2, max(p)):
        q = p[d + 1] / p[d]
        assert q == q_kernel[d]
        product *= q
        assert product == p[d + 1]
        checks.append({"d": d, "qbar": str(q), "required_C": float(q * (d + 1) ** 2)})
    return {"n": n, "A1": str(a1), "checks": checks,
            "max_required_C": max(x["required_C"] for x in checks)}


result = {"status": "passed", "scope": "finite exact cross-check only",
          "results": [one(n) for n in (3, 4, 5)]}
print(json.dumps(result, indent=2))
