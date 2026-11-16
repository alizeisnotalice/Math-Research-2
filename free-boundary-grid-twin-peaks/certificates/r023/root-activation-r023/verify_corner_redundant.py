#!/usr/bin/env python3
"""Finite checks for the corner-redundant subclass and C_3 pattern."""
import json
from fractions import Fraction
from pathlib import Path

B = Path(__file__).resolve().parent
for n in (3, 4, 5):
    a = json.loads((B / f"corner-n{n}.json").read_text())
    c3 = next(x for x in a["by_C_size"] if x["C_size"] == 3)["minimum_p_act"]
    assert Fraction(c3["p_act"]) == Fraction(2, n * n - 4)
    terms = c3["neighbour_completion_terms"]
    assert len(terms) == 2 and terms[0]["h_S_plus_z"] == terms[1]["h_S_plus_z"]
    assert int(c3["h_S"]) == (n * n - 4) * int(terms[0]["h_S_plus_z"])
    full = json.loads((B / f"n{n}-cpp.json").read_text())
    for r, s in zip(a["by_C_size"], full["by_C_size"]):
        k = r["C_size"]
        if k >= 3:
            target = Fraction(2, n * n - k - 1)
            assert Fraction(r["minimum_p_act"]["p_act"]) == target
            assert Fraction(s["minimum_p_act"]["p_act"]) == target

for n in (4, 5):
    restricted = json.loads((B / f"corner-n{n}.json").read_text())
    full = json.loads((B / f"n{n}-cpp.json").read_text())
    assert [x["C_size"] for x in restricted["by_C_size"]] == [x["C_size"] for x in full["by_C_size"]]
    for a, b in zip(restricted["by_C_size"], full["by_C_size"]):
        for key in ("minimum_p_act", "minimum_p_times_C_plus_1", "minimum_p_times_sqrt_C"):
            assert Fraction(a[key]["p_act"]) == Fraction(b[key]["p_act"])
print("Corner C_3 pattern and n4/n5 all-layer restricted/full extremum agreement verified.")
