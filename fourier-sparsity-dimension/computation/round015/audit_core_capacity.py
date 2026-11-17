#!/usr/bin/env python3
"""Exact small-domain audit of kappa(T) <= kappa(C) + kappa(B).

Inputs are the already completed R14 task/result JSONL files.  No lifted-code
capacity DP is invoked here: for each distinct core B this file computes its
support-partition capacity directly on its original <=7 point domain using
integer/F2 operations, then joins it back to every completed T record.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def popcount(value: int) -> int:
    """Python 3.9-compatible exact population count."""
    return bin(value).count("1")


def span(basis: tuple[int, ...]) -> list[int]:
    words = [0]
    for g in basis:
        words += [x ^ g for x in words]
    if len(set(words)) != len(words):
        raise ValueError("dependent core basis")
    return words


def capacity(spec: tuple[int, tuple[int, ...]]) -> dict:
    """Integer recurrence used by the R13 DP, on a small original domain."""
    points, basis = spec
    words = span(basis)
    full = (1 << points) - 1
    try:
        root = words.index(full)
    except ValueError as exc:
        raise ValueError("core B is not unital") from exc
    m = len(words)
    sub_dim = [0] * m
    for s, support in enumerate(words):
        count = sum((cell & ~support) == 0 for cell in words)
        if count == 0 or count & (count - 1):
            raise ValueError("nonlinear containment count")
        sub_dim[s] = count.bit_length() - 1
    order = sorted(range(m), key=lambda s: (popcount(words[s]), s))
    dp = [-10**9] * m
    pick = [0] * m
    dp[0] = 0
    for s in order:
        if s == 0:
            continue
        support = words[s]
        anchor = support & -support
        best = -10**9
        for a, cell in enumerate(words):
            if a and (cell & ~support) == 0 and cell & anchor:
                value = sub_dim[a] + 1 + dp[s ^ a]
                if value > best:
                    best, pick[s] = value, a
        if best < 0:
            raise ValueError("recurrence failed")
        dp[s] = best
    partition: list[int] = []
    remaining = root
    union = 0
    while remaining:
        a = pick[remaining]
        if a == 0 or union & words[a]:
            raise ValueError("invalid backtrace")
        union |= words[a]
        partition.append(words[a])
        remaining ^= a
    if union != full:
        raise ValueError("partition does not cover full support")
    check = sum(
        (sum((h & ~cell) == 0 for h in words).bit_length())
        for cell in partition
    )
    if check != dp[root]:
        raise ValueError("partition capacity check failed")
    return {"points": points, "basis": list(basis), "dimension": len(basis),
            "kappa_b": dp[root], "partition": [hex(cell) for cell in partition]}


def read_jsonl(path: Path) -> list[dict]:
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=15)
    args = ap.parse_args()
    if not 1 <= args.workers <= 15:
        raise SystemExit("workers must be 1..15")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    tasks = read_jsonl(args.tasks)
    task_by_id = {x["id"]: x for x in tasks}
    if len(task_by_id) != len(tasks):
        raise ValueError("duplicate task ID")
    results = read_jsonl(args.results)
    result_by_id = {x["id"]: x for x in results}
    if len(result_by_id) != len(results) or set(result_by_id) != set(task_by_id):
        raise ValueError("task/result coverage mismatch")
    for i, result in result_by_id.items():
        task = task_by_id[i]
        if result["task"] != task or result["returncode"] != 0 or result["status"] != "ok":
            raise ValueError(f"invalid completed T result {i}")
    unique_specs = sorted({(x["points"], tuple(x["B_basis"])) for x in tasks})
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        core_records = list(pool.map(capacity, unique_specs, chunksize=16))
    capacity_by_spec = {(x["points"], tuple(x["basis"])): x for x in core_records}
    comparison_path = args.out_dir / "per-task-comparison.jsonl"
    equality = strict = violations = 0
    gap_hist: dict[str, int] = {}
    smallest: dict | None = None
    with comparison_path.open("w") as handle:
        for i in sorted(task_by_id):
            task, result = task_by_id[i], result_by_id[i]
            core = capacity_by_spec[(task["points"], tuple(task["B_basis"]))]
            rhs = task["base_capacity"] + core["kappa_b"]
            gap = rhs - result["kappa"]
            row = {"id": i, "points": task["points"],
                   "lift_dimension": task["lift_dimension"],
                   "kappa_t": result["kappa"], "kappa_c": task["base_capacity"],
                   "kappa_b": core["kappa_b"], "rhs": rhs, "gap": gap,
                   "equality": gap == 0}
            handle.write(json.dumps(row, separators=(",", ":")) + "\n")
            gap_hist[str(gap)] = gap_hist.get(str(gap), 0) + 1
            if gap == 0:
                equality += 1
            elif gap > 0:
                strict += 1
            else:
                violations += 1
                candidate = {"task": task, "completed_T_result": result,
                             "exact_core_B": core, "comparison": row}
                if smallest is None or (task["points"], task["lift_dimension"], i) < (
                    smallest["task"]["points"], smallest["task"]["lift_dimension"], smallest["task"]["id"]):
                    smallest = candidate
    counter_path = args.out_dir / "minimal-counterexample.json"
    counter_path.write_text(json.dumps({"exists": smallest is not None,
        "ordering": "(points,lift_dimension,id)", "counterexample": smallest},
        indent=2, sort_keys=True) + "\n")
    summary = {
        "schema": "a28-r015-core-capacity-audit-v1",
        "status": "passed" if violations == 0 else "counterexample_found",
        "started_utc": started_utc,
        "tasks": len(tasks), "unique_core_B_codes": len(unique_specs),
        "workers": args.workers, "equality": equality, "strict": strict,
        "violations": violations, "gap_histogram": gap_hist,
        "minimal_counterexample_path": str(counter_path),
        "source_sha256": sha256(Path(__file__)),
        "tasks_sha256": sha256(args.tasks), "results_sha256": sha256(args.results),
        "comparison_sha256": sha256(comparison_path),
        "counterexample_sha256": sha256(counter_path),
        "python": sys.version, "platform": platform.platform(),
        "elapsed_seconds": time.monotonic() - started,
        "scope": "Exact capacity recurrence for all supplied R14 core B codes and comparison with already completed T capacities; no H=C+wB scan.",
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ["status", "tasks", "unique_core_B_codes", "equality", "strict", "violations", "elapsed_seconds"]}, sort_keys=True))
    return 0 if violations == 0 else 10


if __name__ == "__main__":
    raise SystemExit(main())
