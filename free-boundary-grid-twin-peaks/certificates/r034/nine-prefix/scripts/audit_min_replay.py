#!/usr/bin/env python3
"""Check the finite gate and 16-host-edge min replay artifacts."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--work", type=Path, required=True)
ap.add_argument("--source-sha", required=True)
ap.add_argument("--json-out", type=Path)
ns = ap.parse_args()
for mode in ("gate-n7-r45", "gate-n9-direct"):
    row = json.loads((ns.work / "gates" / f"{mode}.json").read_text())
    assert row["status"] == "passed" and row["source_sha256"] == ns.source_sha
row = json.loads((ns.work / "results" / "shard-l1-min-l1-0-16.json").read_text())
assert row["status"] == "completed" and row["graphs_scanned"] == 16 and row["max_abs_R27"] == 0 and row["source_sha256"] == ns.source_sha
result = {"status": "passed", "replayed_modes": ["gate-n7-r45", "gate-n9-direct", "shard-r27-left1-0-16"]}
text = json.dumps(result, indent=2, sort_keys=True) + "\n"
if ns.json_out:
    ns.json_out.write_text(text)
print(text, end="")
