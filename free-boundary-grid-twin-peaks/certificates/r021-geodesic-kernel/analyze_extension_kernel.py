#!/usr/bin/env python3
"""Turn integer extension-kernel output into exact ratios and required local C."""
import json
import sys
from fractions import Fraction
from pathlib import Path

raw = json.loads(Path(sys.argv[1]).read_text())
SOURCE = Path(__file__).resolve().parent / "reference"
reference_path = SOURCE / f"geodesic-cpp-n{raw['n']}.json"
ref_p = None
if reference_path.exists():
    reference = json.loads(reference_path.read_text())
    ref_p = {int(r["d"]): Fraction(int(r["numerator"]), int(reference["A1"]))
             for r in reference["rows"]}
else:
    old = json.loads((SOURCE / "geodesic-r001.json").read_text())
    reference = next(x for x in old["results"] if x["n"] == raw["n"])
    ref_p = {int(r["d"]): Fraction(int(r["numerator"]), int(reference["A1"]))
             for r in reference["rows"]}
local_rows = {}
for line in Path(sys.argv[2]).read_text().splitlines():
    d_text, path, num, den = line.split("\t")
    d = int(d_text)
    q = Fraction(int(num), int(den))
    old = local_rows.get(d)
    if old is None or q > old[0]:
        local_rows[d] = (q, path, num, den)
for row in raw["rows"]:
    j = int(row["prefix_vertices"])
    q = Fraction(int(row["extension_completion_weight"]), int(row["prefix_completion_weight"]))
    if ref_p is not None:
        assert q == ref_p[j + 1] / ref_p[j]
    row["geodesic_extension_probability"] = str(q)
    row["geodesic_extension_decimal"] = float(q)
    row["required_C_for_q_le_C_over_(j+1)^2"] = float(q * (j+1)**2)
    mq, path, num, den = local_rows[j]
    row["max_local_extension_probability"] = str(mq)
    row["max_local_extension_decimal"] = float(mq)
    row["max_local_required_C"] = float(mq * (j+1)**2)
    row["max_local_witness"] = [int(x) for x in path.split(",")]
raw["crosscheck"] = "all weighted extension ratios equal certified p_(j+1)/p_j"
print(json.dumps(raw, indent=2))
