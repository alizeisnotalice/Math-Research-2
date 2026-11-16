#!/usr/bin/env python3
"""Exact consistency checks for first_merge_spatial_dp JSON outputs."""
import json
from fractions import Fraction
from pathlib import Path

EXPECTED = {
    3: (143808, 633600),
    4: (779299174080, 5272309842240),
    5: (7063261481362958280000, 67631454448823894520000),
}

rows = []
for n, (a2, distance_sum) in EXPECTED.items():
    data = json.loads(Path(f"first-merge-n{n}.json").read_text())
    assert int(data["A2_from_W"]) == a2
    assert int(data["D_euclidean_squared_from_W"]) == distance_sum
    ledger = data["rank_minbranch_span"]
    assert sum(int(x["W"]) for x in ledger) == a2
    assert sum(int(x["D_euclidean_squared"]) for x in ledger) == distance_sum
    assert all(x["s"] >= 3 and x["min_branch"] >= 1
               and x["branch_span_manhattan"] >= 2 for x in ledger)
    rows.append({
        "n": n,
        "A2": a2,
        "D_euclidean_squared": distance_sum,
        "Q": str(Fraction(distance_sum, n*n*a2)),
        "strata": len(ledger),
    })
Path("verify-first-merge-spatial.json").write_text(
    json.dumps({"status": "passed", "rows": rows}, indent=2) + "\n"
)
print(json.dumps({"status": "passed", "rows": rows}, indent=2))
