#!/usr/bin/env python3
"""Summarize forced-rank CPU SIS output; diagnostics only."""
import json
import math
import sys
from pathlib import Path

out = Path(sys.argv[1])
rows = [json.loads(x) for x in (out / "summary.jsonl").read_text().splitlines()]
n = rows[0]["n"]
def f(r, key):
    return float(r[key])
a2 = sum(f(r, "A2_mean") for r in rows)
e2 = sum(f(r, "A2_mean") * f(r, "Q_e2") * n*n for r in rows)
rho2 = sum(f(r, "A2_mean") * f(r, "rho2_mean") for r in rows)
half = sum(f(r, "A2_mean") * f(r, "far_half_weight") for r in rows)
quarter = sum(f(r, "A2_mean") * f(r, "far_quarter_weight") for r in rows)
result = {"n": n, "strata": len(rows), "A2_sum": a2,
          "Q_e2_ratio": e2/(n*n*a2) if a2 else None,
          "rho2_ratio": rho2/a2 if a2 else None,
          "far_half_ratio": half/a2 if a2 else None,
          "far_quarter_ratio": quarter/a2 if a2 else None,
          "A1_from_tau2": f(rows[0], "A1_mean"),
          "min_A2_ess": min(f(r,"A2_ess") for r in rows if f(r,"A2_mean")),
          "max_A2_top_share": max(f(r,"A2_top_share") for r in rows if f(r,"A2_mean")),
          "dead_total": sum(r["dead"] for r in rows),
          "samples_per_stratum": rows[0]["samples"],
          "numerical_only": True}
(out / "aggregate.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
