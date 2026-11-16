#!/usr/bin/env python3
import csv
import glob
import json
import math
import pathlib

here = pathlib.Path(__file__).resolve().parent

checks = []
for old_path in sorted(here.glob("old-n*.csv")):
    new_path = here / old_path.name.replace("old-", "new-", 1)
    old_rows = list(csv.DictReader(old_path.open()))
    new_rows = list(csv.DictReader(new_path.open()))
    n = int(old_rows[0]["n"])
    factorial = math.factorial(n * n - 1)
    same = len(old_rows) == len(new_rows)
    max_relative_normalizer_error = 0.0
    for old, new in zip(old_rows, new_rows):
        comparable = ("seed", "min_parent_ess", "dead", "accept", "proposals")
        same &= int(old["seed"]) == int(new["seed"])
        same &= int(old["dead"]) == int(new["dead"])
        for field in ("min_parent_ess", "accept", "proposals"):
            same &= float(old[field]) == float(new[field])
        lhs = float(old["Z"])
        rhs = float(new["probZ"]) * factorial
        err = abs(lhs-rhs) / max(1.0, abs(lhs))
        max_relative_normalizer_error = max(max_relative_normalizer_error, err)
        same &= err < 2e-14
    checks.append({"old": old_path.name, "new": new_path.name,
                   "rows": len(old_rows), "all_shared_fields_identical": same,
                   "max_relative_normalizer_error": max_relative_normalizer_error})

terminal = []
for path in sorted(here.glob("n16-p*.csv")):
    rows = list(csv.DictReader(path.open()))
    if not rows:
        continue
    terminal.append({"file": path.name, "rows": len(rows),
        "all_nonzero": all(float(x["probZ"]) > 0 for x in rows),
        "all_zero_dead": all(int(x["dead"]) == 0 for x in rows),
        "min_ess": min(float(x["min_parent_ess"]) for x in rows),
        "proposal_total": sum(int(x["proposals"]) for x in rows),
        "accept_total": sum(int(x["accept"]) for x in rows),
        "terminal_hashes": [x["terminal_hash"] for x in rows]})

report = {"old_new_comparison": checks, "n16_terminal_assertion_runs": terminal}
print(json.dumps(report, indent=2, sort_keys=True))
if not all(x["all_shared_fields_identical"] for x in checks):
    raise SystemExit("old/new fixed-root comparisons failed")
if terminal and not all(x["all_nonzero"] and x["all_zero_dead"] for x in terminal):
    raise SystemExit("n16 terminal checks failed")
