#!/usr/bin/env python3
import csv
import json
import math
import pathlib
import statistics

here = pathlib.Path(__file__).resolve().parent
path = here / "n16-root119-p65536-r64-l16-direct.csv"
rows = list(csv.DictReader(path.open()))
assert len(rows) == 64
assert sorted(int(x["rep"]) for x in rows) == list(range(64))
assert len({int(x["seed"]) for x in rows}) == 64
assert all(int(x["n"]) == 16 and int(x["N"]) == 256 and int(x["root"]) == 119 for x in rows)
assert all(int(x["L"]) == 16 and int(x["dead"]) == 0 for x in rows)
base = 0xA0249000
for row in rows:
    expected = base + 119 * 1000000007 + int(row["rep"]) * 1009 + 16 * 1000003
    assert int(row["seed"]) == expected
z = [float(x["probZ"]) for x in rows]
assert all(math.isfinite(x) and x > 0 for x in z)
assert all(math.isfinite(float(x["logZ"])) for x in rows)
mean = statistics.fmean(z)
se = statistics.stdev(z) / math.sqrt(len(z))
out = {
    "kind": "finite_SMC_diagnostic_not_certified_count",
    "source": "fixed_root_n16.cpp",
    "parameters": {"n": 16, "root": 119, "root_xy": [7, 7], "P": 65536,
                   "replicates": 64, "move_steps": 16, "threads": 12,
                   "seed_base": "0xA0249000"},
    "strict_checks": {"row_count": len(rows), "replicates_exactly_0_to_63": True,
                      "seeds_unique_and_formula_match": True, "all_dead_zero": True,
                      "finite_positive_probZ": True,
                      "terminal_peak_assertions_were_enabled_in_binary": True},
    "estimate": {"mean_probZ": mean, "sample_standard_error": se,
                 "relative_sample_standard_error": se / mean,
                 "min_probZ": min(z), "max_probZ": max(z),
                 "mean_min_parent_ess": statistics.fmean(float(x["min_parent_ess"]) for x in rows),
                 "min_min_parent_ess": min(float(x["min_parent_ess"]) for x in rows),
                 "mean_acceptance": statistics.fmean(float(x["accept"]) / float(x["proposals"]) for x in rows),
                 "unique_terminal_hashes": len({x["terminal_hash"] for x in rows})},
}
print(json.dumps(out, indent=2, sort_keys=True))
