#!/usr/bin/env python3
"""Validate and compact JSONL output from gpu_original_law_sis.swift."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import re


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    args = ap.parse_args()
    data = [json.loads(line) for line in (args.run_dir / "results.jsonl").read_text().splitlines()]
    assert data and data[0]["kind"] == "metadata"
    groups = defaultdict(list)
    checks = []
    batches = defaultdict(list)
    for row in data[1:]:
        if row["kind"] == "trajectory_check":
            checks.append(row)
        elif row["kind"] == "batch":
            batches[(row["n"], row["mode"])].append(row)
        elif row["kind"] == "summary":
            groups[(row["n"], row["mode"])].append(row)
        else:
            raise AssertionError(row["kind"])
    assert all(x["ok"] is True for x in checks)
    assert all(len(v) == 1 for v in groups.values())
    summary = []
    for key, entries in sorted(groups.items()):
        row = entries[0]
        bs = batches[key]
        assert len(bs) * data[0]["batch"] == row["samples"]
        assert [b["batch_index"] for b in bs] == list(range(len(bs)))
        assert all(b["gpu_seconds"] >= 0 and b["wall_seconds"] >= 0 for b in bs)
        assert row["invalid"] == 0
        summary.append({
            "n": row["n"], "mode": row["mode"], "samples": row["samples"],
            "A_ratio_to_exact": row["A_ratio_to_exact"],
            "Q_euclidean_hat": row["Q_euclidean_hat"],
            "Q_error_to_exact": row["Q_error_to_exact"],
            "ESS_A1": row["ESS_A1"], "ESS_A2": row["ESS_A2"],
            "top_weight_A1": row["top_weight_A1"], "top_weight_A2": row["top_weight_A2"],
            "gpu_seconds": row["gpu_seconds"], "submit_wait_seconds": row["submit_wait_seconds"],
            "batch_count": len(bs),
        })
    agx = args.run_dir / "agx-samples.txt"
    util = []
    if agx.exists():
        util = [int(x) for x in re.findall(r'Device Utilization %"=(\d+)', agx.read_text())]
    out = {"status": "complete", "metadata": data[0], "trajectory_checks": len(checks),
           "summaries": summary,
           "agx_device_samples": len(util),
           "agx_device_mean": sum(util) / len(util) if util else None,
           "agx_device_ge_70": sum(x >= 70 for x in util) if util else None,
           "agx_device_max": max(util) if util else None,
           "limits": "finite-precision deterministic RNG and self-normalized Q; numerical exploration only, not a certified conditional probability estimate"}
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
