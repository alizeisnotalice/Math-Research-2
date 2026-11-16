#!/usr/bin/env python3
"""Exact finite audit of the completion-weighted global geodesic recursion.

Inputs are the certified integer global prefix numerators.  For consecutive
d, qbar_d=p_{d+1}/p_d is the h-weighted average of the outward geodesic
extension kernel, not a maximum over individual prefixes.
"""
from fractions import Fraction
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = Path("/Volumes/data/CodexTemp/T2/a022-main-rounds/free-boundary-grid-twin-peaks/computation/results")
INDEPENDENT_KERNEL = Path("/Volumes/data/CodexTemp/T2/a022-geodesic-average-r020")


def rows_for(n):
    cpp = SOURCE / f"geodesic-cpp-n{n}.json"
    if cpp.exists():
        raw = json.loads(cpp.read_text())
        return int(raw["A1"]), raw["rows"]
    old = json.loads((SOURCE / "geodesic-r001.json").read_text())
    raw = next(x for x in old["results"] if x["n"] == n)
    return int(raw["A1"]), raw["rows"]


def one(n):
    a1, rows = rows_for(n)
    p = {int(r["d"]): Fraction(int(r["numerator"]), a1) for r in rows}
    checks = [{"d": 1, "p_d": str(p[1]), "p_next": str(p[2]),
               "normalization_note": "p_2=1; F_2/F_1=1/2 because singleton paths have one, not two, orientations"}]
    kernel = json.loads((INDEPENDENT_KERNEL / f"n{n}-extension-analysis.json").read_text())
    kernel_q = {int(r["prefix_vertices"]): Fraction(r["geodesic_extension_probability"])
                for r in kernel["rows"]}
    for d in range(2, max(p)):
        q = p[d + 1] / p[d]
        assert q == kernel_q[d]
        checks.append({"d": d, "p_d": str(p[d]), "p_next": str(p[d + 1]),
                       "qbar": str(q), "qbar_decimal": float(q),
                       "required_C_for_qbar": float(q * (d + 1) ** 2)})
    product = Fraction(1)
    for row in checks[1:]:
        product *= Fraction(row["qbar"])
        assert product == p[row["d"] + 1]
    return {"n": n, "A1": str(a1), "rows": checks,
            "independent_kernel_crosscheck": "passed exact qbar equality for every d >= 2",
            "maximum_required_average_C": max(x["required_C_for_qbar"] for x in checks[1:])}


out = {"definition": "qbar_d=p_(d+1)/p_d, the completion-weighted global outward-extension average",
       "warning": "finite exact diagnostics only; neither a uniform proof nor a counterexample",
       "results": [one(n) for n in (3, 4, 5)]}
(HERE / "weighted-geodesic-recursion.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
