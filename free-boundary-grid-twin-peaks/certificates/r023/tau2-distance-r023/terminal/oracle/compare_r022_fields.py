#!/usr/bin/env python3
"""Verify r023 changed no pre-existing replicate field except timings."""
import json
import pathlib
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: compare_r022_fields.py r022.jsonl r023.jsonl")


def read(path):
    return next(json.loads(line) for line in pathlib.Path(path).read_text().splitlines()
                if json.loads(line)["kind"] == "replicate")


old, new = read(sys.argv[1]), read(sys.argv[2])
for key in ("gpu_seconds", "wall_seconds", "tau2_distance_L_hist_groups",
            "tau2_distance_L2_groups", "tau2_distance_scaled_mean_L_hist",
            "tau2_distance_scaled_mean_L2"):
    old.pop(key, None)
    new.pop(key, None)
assert old == new, {"old": old, "new": new}
print(json.dumps({"passed": True, "fields": sorted(old)}))
