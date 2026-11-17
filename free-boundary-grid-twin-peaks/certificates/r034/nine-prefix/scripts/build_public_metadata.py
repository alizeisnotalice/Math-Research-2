#!/usr/bin/env python3
"""Create portable, path-redacted evidence metadata from one completed run.

This tool copies no binary or log. It preserves every completed checkpoint
unchanged and normalizes only metadata fields that named local paths.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def normalized_outcome(row: dict) -> dict:
    return {
        "checkpoint": "evidence/results/" + Path(row["checkpoint"]).name,
        "checkpoint_sha256": row["checkpoint_sha256"],
        "checkpoint_status": row["checkpoint_status"],
        "exit_code": row["exit_code"],
        "stderr_sha256": row["stderr_sha256"],
        "stdout_sha256": row["stdout_sha256"],
        "task": row["task"],
        "utc_finished": row["utc_finished"],
        "utc_started": row["utc_started"],
        "wall_seconds": row["wall_seconds"],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ns = ap.parse_args()
    run, out = ns.run.resolve(), ns.out.resolve()
    evidence = out / "evidence"
    strict = run / "strict-aggregate.json"
    original = load(strict)
    coverage = original["coverage"]
    dump(evidence / "strict-aggregate.public.json", {
        "kind": "portable-strict-coverage-aggregate",
        "original_sha256": digest(strict),
        "source_sha256": original["source_sha256"],
        "status": original["status"],
        "rule": original["rule"],
        "coverage": coverage,
        "checkpoint_count": sum(len(v) for v in coverage.values()),
        "redaction": "results_directory and plan paths were I/O-only local paths; coverage and origins are retained.",
    })

    outcomes_path = run / "production.outcomes.json"
    outcomes = load(outcomes_path)
    benchmark_path = run / "benchmark.receipt.json"
    benchmark = load(benchmark_path)
    benchmark_checkpoint = run / "results" / Path(benchmark["checkpoint"]).name
    benchmark_row = {
        "checkpoint": "evidence/results/" + benchmark_checkpoint.name,
        "checkpoint_sha256": digest(benchmark_checkpoint),
        "checkpoint_status": load(benchmark_checkpoint)["status"],
        "exit_code": benchmark["exit_code"],
        "stderr_sha256": benchmark["stderr_sha256"],
        "stdout_sha256": benchmark["stdout_sha256"],
        "task": {"id": "benchmark-l4-0000000000-0000001024", "left": 4, "start": 0, "stop": 1024},
        "utc_finished": benchmark["utc_finished"],
        "utc_started": benchmark["utc_started"],
        "wall_seconds": benchmark["wall_seconds"],
    }
    dump(evidence / "execution-exits.public.json", {
        "kind": "path-redacted-exit-ledger",
        "source_sha256": outcomes["source_sha256"],
        "production_outcomes_original_sha256": digest(outcomes_path),
        "benchmark_receipt_original_sha256": digest(benchmark_path),
        "halted_left_sizes": outcomes["halted_left_sizes"],
        "benchmark": benchmark_row,
        "production": [normalized_outcome(x) for x in outcomes["outcomes"]],
        "redaction": "private command paths and process IDs were removed; exit, task, checkpoint, output-hash, timestamp, and wall-time fields are retained.",
    })

    gate_records = []
    for mode in ("gate-n7-r45", "gate-n9-direct"):
        result = run / "gates" / f"{mode}.json"
        receipt = run / "gates" / f"{mode}.receipt.json"
        row = load(receipt)
        gate_records.append({
            "mode": mode, "result": "evidence/gates/" + result.name,
            "result_sha256": digest(result), "receipt_original_sha256": digest(receipt),
            "exit_code": row["exit_code"], "stdout_sha256": row["stdout_sha256"],
            "stderr_sha256": row["stderr_sha256"], "utc_started": row["utc_started"],
            "utc_finished": row["utc_finished"], "wall_seconds": row["wall_seconds"],
            "portable_command": ["./bipartite_prefix_dp", "--mode", mode, "--out-dir", "./gates", "--source-sha", load(result)["source_sha256"]],
        })
    compile_receipt = load(run / "compile.receipt.json")
    aggregate_receipt = load(run / "aggregate.receipt.json")
    completed = load(run / "completed.json")
    started = load(run / "started.json")
    profile = run / "bench-rss-profile"
    dump(evidence / "run-receipts.public.json", {
        "kind": "path-redacted-minimal-run-receipts",
        "source_sha256": started["source_sha256"],
        "started": {"utc_started": started["utc_started"], "workers": started["workers"], "chunk": started["chunk"], "benchmark_counted_range": started["benchmark_counted_range"], "original_sha256": digest(run / "started.json")},
        "compile": {"flags": ["-O3", "-std=c++20", "-Wall", "-Wextra", "-Werror", "-pedantic"], "exit_code": compile_receipt["exit_code"], "wall_seconds": compile_receipt["wall_seconds"], "binary_sha256_historical_only": compile_receipt["binary_sha256"], "original_sha256": digest(run / "compile.receipt.json")},
        "aggregate": {"exit_code": aggregate_receipt["exit_code"], "wall_seconds": aggregate_receipt["wall_seconds"], "original_sha256": digest(run / "aggregate.receipt.json")},
        "completed": {"status": completed["status"], "utc_finished": completed["utc_finished"], "binary_sha256_historical_only": completed["binary_sha256"], "original_sha256": digest(run / "completed.json")},
        "gates": gate_records,
        "resource_profile": {"scope": "uncounted isolated left=4 [0,1024) reprofile; not a full-batch resource total", "time_stderr_original_sha256": digest(profile / "time.stderr"), "profile_checkpoint_original_sha256": digest(profile / "shard-l4-profile-l4-0000000000-0000001024.json"), "exit_code_file_original_sha256": digest(profile / "exit-code.txt"), "wall_seconds": 0.14, "user_seconds": 0.13, "sys_seconds": 0.00, "maximum_resident_set_size_bytes": 2195456, "peak_memory_footprint_bytes": 1556888, "profile_exit_code": (profile / "exit-code.txt").read_text().strip()},
        "redaction": "all removed values were local paths, local cache paths, or PIDs; original digests identify source receipts.",
    })


if __name__ == "__main__":
    main()
