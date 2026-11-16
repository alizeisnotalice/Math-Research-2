#!/usr/bin/env python3
"""Recompute r014 diagnostics from the raw JSONL; no iid claim is made."""
from collections import defaultdict
from pathlib import Path
import json
import math


def summarize(path: Path) -> None:
    groups = defaultdict(list)
    audits = []
    for line in path.read_text().splitlines():
        row = json.loads(line)
        if row["kind"] == "replicate":
            groups[row["shape"]].append(row)
        else:
            audits.append(row)
    assert all(row["kind"] == "cpu_gpu_small_audit" and row["ok"] for row in audits)
    for shape, rows in sorted(groups.items()):
        top = max(row["logZ"] for row in rows)
        weights = [math.exp(row["logZ"] - top) for row in rows]
        total = sum(weights)
        p = sum(weight * row["p"] for weight, row in zip(weights, rows)) / total
        zess = total * total / sum(weight * weight for weight in weights)
        result = {
            "shape": shape,
            "replicates": len(rows),
            "particles": rows[0]["particles"],
            "scaled_Z_weighted": rows[0]["boundary"] * p,
            "exact": rows[0]["exact_scaled"],
            "Z_ESS": zess,
            "top_Z_share": max(weights) / total,
            "min_roots": min(row["min_roots"] for row in rows),
            "metal_gpu_seconds": sum(row["gpu_seconds"] for row in rows),
            "submit_wait_seconds": sum(row["wall_seconds"] for row in rows),
        }
        print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    summarize(Path(__file__).with_name("production.jsonl"))
