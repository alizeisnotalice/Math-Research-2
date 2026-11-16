#!/usr/bin/env python3
"""Assert all pre-r022 replicate fields agree exactly, except timing."""
import json
import pathlib
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: compare_legacy_replicate.py old.jsonl new.jsonl")

def replicate(path):
    return next(json.loads(line) for line in pathlib.Path(path).read_text().splitlines()
                if json.loads(line)["kind"] == "replicate")

old, new = replicate(sys.argv[1]), replicate(sys.argv[2])
for key in ("gpu_seconds", "wall_seconds",
            "geodesic_prefix_group_counts", "scaled_geodesic_prefix_p"):
    old.pop(key, None)
    new.pop(key, None)
assert old == new, {"old": old, "new": new}
print(json.dumps({"passed": True, "compared_fields": sorted(old)}))
