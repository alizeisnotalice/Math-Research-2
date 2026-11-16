#!/usr/bin/env python3
"""Compare the independent C++ n3/n4 extrema with the Python exact scan."""
import json
from fractions import Fraction
from pathlib import Path

BASE = Path(__file__).resolve().parent
old = json.loads((BASE / "root-activation-scan-n3-n4.json").read_text())["results"]
for n in (3, 4):
    got = json.loads((BASE / f"n{n}-cpp.json").read_text())
    want = old[n - 3]
    assert got["eligible_pairs"] == want["eligible_pairs"]
    for x, y in zip(got["by_C_size"], want["by_C_size"]):
        assert x["C_size"] == y["C_size"]
        for key in ("minimum_p_act", "minimum_p_times_C_plus_1", "minimum_p_times_sqrt_C"):
            assert Fraction(x[key]["p_act"]) == Fraction(y[key]["p_act"])
            if key == "minimum_p_times_C_plus_1":
                assert Fraction(x[key]["p_times_C_plus_1"]) == Fraction(y[key]["p_times_C_plus_1"])
print("C++ and Python n3/n4 extrema agree at every C-size (D4 tie representatives may differ).")
