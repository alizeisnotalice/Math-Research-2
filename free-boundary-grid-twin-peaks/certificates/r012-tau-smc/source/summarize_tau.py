#!/usr/bin/env python3
"""Aggregate tau strata without converting ratios into a proof."""
from collections import defaultdict
from pathlib import Path
import json
import math
import sys


rows = defaultdict(list)
for line in Path(sys.argv[1]).read_text().splitlines():
    x = json.loads(line)
    if x.get("kind") == "replicate":
        rows[x["tau"]].append(x)

sum_z = 0.0
sum_u = 0.0
for tau, xs in sorted(rows.items()):
    z = sum(x["Z"] for x in xs) / len(xs)
    u = sum(x["U"] for x in xs) / len(xs)
    sum_z += z
    sum_u += u
    logz = [x["logZ"] for x in xs if x["logZ"] is not None and math.isfinite(x["logZ"])]
    w = [math.exp(v - max(logz)) for v in logz] if logz else []
    ess = sum(w) ** 2 / sum(v*v for v in w) if w else 0.0
    print(json.dumps({"tau": tau, "reps": len(xs), "Z_mean": z, "U_mean": u,
                      "dead": sum(x["dead"] for x in xs),
                      # Older files used the misleading key min_distinct_ancestors.
                      # In both formats this is only a direct-parent diagnostic.
                      "min_one_step_distinct_parents": min(x.get("min_one_step_distinct_parents", x.get("min_distinct_ancestors")) for x in xs),
                      "replicate_Z_ESS": ess}, sort_keys=True))
print(json.dumps({"sum_Z_mean": sum_z, "sum_U_mean": sum_u,
                  "ratio_numerical_only": sum_u / sum_z if sum_z else None}, sort_keys=True))
