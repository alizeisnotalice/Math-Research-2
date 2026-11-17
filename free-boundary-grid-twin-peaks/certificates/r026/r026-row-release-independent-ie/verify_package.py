#!/usr/bin/env python3
"""Verify the frozen selected-prefix IE outputs and the second local-minima IE."""
import json
from pathlib import Path

root = Path(__file__).parent
expected = {
    3: ("156", "6"),
    4: ("4213440", "5600"),
    5: ("153450311313600", "3971406600"),
    6: ("20779014929894865426585600", "6637177983794054400"),
}
checked = []
for n, (ht, hr) in expected.items():
    t = json.loads((root / f"n{n}-T.json").read_text())
    r = json.loads((root / f"n{n}-R.json").read_text())
    assert t["hS"] == ht and r["hS"] == hr
    if n <= 5:
        assert t["cross_hS_equal"] and r["cross_hS_equal"]
    checked.append({"n": n, "H_T": ht, "H_R": hr,
                    "direct_cross": n <= 5})
alternative = {x["case"]: x["value"] for x in json.loads((root / "run-n6.json").read_text())["rows"]}
assert alternative == {"HT": expected[6][0], "HR": expected[6][1]}
print(json.dumps({"all_match": True, "selected_prefix_ie": checked,
                  "second_algorithm_n6_match": True}, indent=2))
