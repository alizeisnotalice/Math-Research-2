#!/usr/bin/env python3
"""Compare the n=7 inverse-fibre SMC aggregate with existing exact totals."""
import csv
from fractions import Fraction
import json
from pathlib import Path


BASE = Path(__file__).parent
ROOT = Path("/Volumes/data/CodexTemp/T2/a022-main-rounds/free-boundary-grid-twin-peaks/"
            "certificates/r018-root-moments/n7-d4-root-moments.json")
with (BASE / "n7-p4096-r32-l16.csv").open() as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 32
assert all(int(r["tau"]) == 50 and int(r["dead"]) == 0 for r in rows)
z = sum(float(r["Z"]) for r in rows)
zm2 = sum(float(r["Zmean_m2"]) for r in rows)
zm3 = sum(float(r["Zmean_m3"]) for r in rows)
mean_z, ratio2, ratio3 = z/len(rows), zm2/z, zm3/z
root = json.loads(ROOT.read_text())
a1, d2, d3 = int(root["A1"]), int(root["distance_counts"]["2"]), int(root["distance_counts"]["3"])
exact2, exact3 = Fraction(d2, a1), Fraction(d3, a1)
errors = [abs(mean_z/a1-1), abs(ratio2/float(exact2)-1), abs(ratio3/float(exact3)-1)]
assert max(errors) < .10, errors
out = {"status": "passed", "scope": "deterministic finite-MC diagnostic gate only",
       "mean_Z": mean_z, "weighted_ratio_m2": ratio2, "weighted_ratio_m3": ratio3,
       "exact_A1": str(a1), "exact_A2d2_over_A1": str(exact2),
       "exact_A2d3_over_A1": str(exact3), "relative_errors": errors,
       "tolerance": .10}
(BASE / "n7-gate-verified.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
