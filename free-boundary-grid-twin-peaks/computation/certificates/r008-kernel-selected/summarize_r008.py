#!/usr/bin/env python3
"""Create a compact exact-integer summary of r008 saved outputs."""
from fractions import Fraction
import glob
import json
from pathlib import Path

base = Path(__file__).parent


def ratio(x):
    return Fraction(int(x["boundary_times_p_numerator"]), int(x["p_denominator"]))


n6 = []
for p in sorted(glob.glob(str(base / "n6-*.json"))):
    x = json.loads(Path(p).read_text())
    n6.append({"file": Path(p).name, "ratio": [str(ratio(x).numerator), str(ratio(x).denominator)],
               "holds": x["candidate_holds"]})

root = base / "n7-critical-d4"
m = json.loads((root / "manifest.json").read_text())
n7 = []
for r in m["runs"]:
    x = json.loads((root / r["json"]).read_text())
    n7.append({"mask": r["mask"], "ratio": [str(ratio(x).numerator), str(ratio(x).denominator)],
               "holds": x["candidate_holds"], "independent_sets": x["independent_sets"]})

n8 = json.loads((base / "n8-extended-L.json").read_text())
out = {
    "n4_cross": json.loads((base / "n4-critical-cross.json").read_text()),
    "n5_cross": json.loads((base / "n5-critical-cross.json").read_text()),
    "n6_selected": n6,
    "n7_d4_critical": {"runs": len(n7), "all_hold": all(x["holds"] for x in n7),
                         "maximum": max(n7, key=lambda x: Fraction(*map(int, x["ratio"]))),
                         "entries": n7},
    "n7_extended_L": json.loads((base / "n7-extended-L-cross.json").read_text()),
    "n8_extended_L": {"output": n8, "strict_counterexample": not n8["candidate_holds"],
                        "positive_gap": str(14 * int(n8["p_numerator"]) - 4 * int(n8["hS"]))}
}
(base / "r008-summary.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
