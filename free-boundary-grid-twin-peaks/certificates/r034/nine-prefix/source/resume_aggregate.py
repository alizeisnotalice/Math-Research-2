#!/usr/bin/env python3
"""Prepare hash-chained suffix resumes and strictly aggregate remote receipts.

This program deliberately never launches a worker.  It recognizes an old
``running``/``failed`` checkpoint only as a certified scanned prefix, never as
completion of its advertised whole range.  A suffix worker must receive a new
shard id because the C++ verifier refuses checkpoint overwrites.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path
import json
from typing import Any

TOTALS = {1: 256, 2: 16384, 3: 262144, 4: 1048576}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict[str, Any]:
    try:
        obj = json.loads(path.read_text())
    except Exception as exc:  # JSON failures must never be treated as zero work.
        raise ValueError(f"unreadable JSON: {path}: {exc}") from exc
    if not isinstance(obj, dict):
        raise ValueError(f"checkpoint is not an object: {path}")
    return obj


def require_source(row: dict[str, Any], source_sha: str, path: Path) -> None:
    if row.get("source_sha256") != source_sha:
        raise ValueError(f"source SHA mismatch: {path}")


def integer(row: dict[str, Any], key: str, path: Path) -> int:
    value = row.get(key)
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"missing/noninteger {key}: {path}")
    return value


def validate_interval(row: dict[str, Any], path: Path, source_sha: str, *, prefix: bool) -> tuple[int, int, int]:
    require_source(row, source_sha, path)
    left = integer(row, "left_size", path)
    start, stop, next_mask = (integer(row, k, path) for k in ("range_start", "range_stop_exclusive", "next_mask"))
    scanned, max_abs = integer(row, "graphs_scanned", path), integer(row, "max_abs_R27", path)
    if left not in TOTALS or not (0 <= start <= next_mask <= stop <= TOTALS[left]):
        raise ValueError(f"invalid interval bounds: {path}")
    expected = (next_mask - start) if prefix else (stop - start)
    if scanned != expected or max_abs != 0 or "counterexample" in row:
        raise ValueError(f"uncertified scanned interval: {path}")
    if not prefix and next_mask != stop:
        raise ValueError(f"completed checkpoint lacks complete interval: {path}")
    return left, start, next_mask if prefix else stop


def atomic_write_new(path: Path, obj: dict[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite plan/aggregate: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    if tmp.exists():
        raise FileExistsError(f"stale temp must be inspected, not replaced: {tmp}")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def cmd_prepare(ns: argparse.Namespace) -> None:
    old = ns.old_checkpoint.resolve()
    row = load(old)
    if row.get("status") not in {"running", "failed"}:
        raise ValueError("only an old running/failed checkpoint can certify a prefix")
    left, start, next_mask = validate_interval(row, old, ns.source_sha, prefix=True)
    stop = integer(row, "range_stop_exclusive", old)
    if next_mask == stop:
        raise ValueError("old checkpoint already covers its entire interval; do not manufacture a resume")
    shard_id = ns.new_shard_id
    if not shard_id or any(ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for ch in shard_id):
        raise ValueError("new shard id must be nonempty ASCII alphanumeric/_/-")
    expected_file = f"shard-l{left}-{shard_id}.json"
    plan = {
        "kind": "hash-chained-prefix-resume-plan",
        "source_sha256": ns.source_sha,
        "old_checkpoint": str(old),
        "old_checkpoint_sha256": digest(old),
        "certified_prefix": {"left_size": left, "start": start, "stop_exclusive": next_mask,
                             "graphs_scanned": next_mask - start, "status": row["status"]},
        "required_new_part": {"left_size": left, "start": next_mask, "stop_exclusive": stop,
                              "new_shard_id": shard_id, "expected_checkpoint_basename": expected_file},
        "remote_command": ["nice", "-n", "10", "./bipartite_prefix_dp", "--mode", "shard-r27",
                           "--left-size", str(left), "--start", str(next_mask), "--stop", str(stop),
                           "--shard-id", shard_id, "--out-dir", ns.out_dir, "--source-sha", ns.source_sha,
                           "--checkpoint-every", str(ns.checkpoint_every)],
        "rule": "Only the recorded [start,next_mask) prefix is certified from the old file. The suffix is zero evidence until its distinct checkpoint is status=completed.",
    }
    atomic_write_new(ns.plan.resolve(), plan)


def plan_prefixes(paths: list[Path], source_sha: str) -> tuple[list[tuple[int, int, int, str]], set[Path]]:
    intervals: list[tuple[int, int, int, str]] = []
    approved_old: set[Path] = set()
    for plan_path in paths:
        plan = load(plan_path)
        if plan.get("kind") != "hash-chained-prefix-resume-plan" or plan.get("source_sha256") != source_sha:
            raise ValueError(f"invalid or source-mismatched plan: {plan_path}")
        old = Path(plan.get("old_checkpoint", "")).resolve()
        if not old.is_file() or digest(old) != plan.get("old_checkpoint_sha256"):
            raise ValueError(f"old checkpoint hash chain broken: {plan_path}")
        row = load(old)
        if row.get("status") not in {"running", "failed"}:
            raise ValueError(f"old status changed after plan: {old}")
        left, start, end = validate_interval(row, old, source_sha, prefix=True)
        cert = plan.get("certified_prefix", {})
        if cert != {"left_size": left, "start": start, "stop_exclusive": end,
                    "graphs_scanned": end-start, "status": row["status"]}:
            raise ValueError(f"plan prefix differs from old receipt: {plan_path}")
        intervals.append((left, start, end, f"old-hash:{digest(old)}"))
        approved_old.add(old)
    return intervals, approved_old


def cmd_aggregate(ns: argparse.Namespace) -> None:
    results = ns.results.resolve()
    plans = [p.resolve() for p in ns.plan]
    intervals, approved_old = plan_prefixes(plans, ns.source_sha)
    if not results.is_dir():
        raise ValueError(f"results directory does not exist: {results}")
    seen_files: list[str] = []
    for path in sorted(results.glob("shard-*.json")):
        row = load(path); status = row.get("status")
        seen_files.append(str(path))
        if status == "counterexample":
            raise ValueError(f"counterexample blocks aggregation: {path}")
        if status == "completed":
            left, start, end = validate_interval(row, path, ns.source_sha, prefix=False)
            intervals.append((left, start, end, f"completed:{digest(path)}"))
        elif status in {"running", "failed"}:
            # The old file is recognized only through its explicit hash-chained plan.
            if path.resolve() not in approved_old:
                raise ValueError(f"unapproved incomplete checkpoint: {path}")
        else:
            raise ValueError(f"unknown shard status: {path}")
    if not seen_files:
        raise ValueError("no shard JSON files found")
    coverage: dict[str, list[dict[str, Any]]] = {}
    for left, total in TOTALS.items():
        rows = sorted((x for x in intervals if x[0] == left), key=lambda x: (x[1], x[2], x[3]))
        cursor = 0
        normalized=[]
        for _, start, end, origin in rows:
            if start != cursor:
                raise ValueError(f"left {left}: gap or overlap at {cursor}; next interval is [{start},{end}) from {origin}")
            if end <= start:
                raise ValueError(f"left {left}: empty/nonpositive interval")
            normalized.append({"start": start, "stop_exclusive": end, "origin": origin})
            cursor = end
        if cursor != total:
            raise ValueError(f"left {left}: incomplete coverage [0,{cursor}) != [0,{total})")
        coverage[str(left)] = normalized
    aggregate = {"status": "completed-strict-coverage", "source_sha256": ns.source_sha,
                 "results_directory": str(results), "plans": [str(x) for x in plans],
                 "coverage": coverage,
                 "rule": "Only completed shard ranges and hash-chained old certified prefixes are counted; no advertised but unscanned suffix is counted."}
    atomic_write_new(ns.out.resolve(), aggregate)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("prepare-resume")
    a.add_argument("--old-checkpoint", type=Path, required=True)
    a.add_argument("--new-shard-id", required=True)
    a.add_argument("--out-dir", required=True, help="remote results directory used in the emitted command")
    a.add_argument("--source-sha", required=True)
    a.add_argument("--checkpoint-every", type=int, default=1000)
    a.add_argument("--plan", type=Path, required=True)
    a.set_defaults(func=cmd_prepare)
    b = sub.add_parser("aggregate")
    b.add_argument("--results", type=Path, required=True)
    b.add_argument("--source-sha", required=True)
    b.add_argument("--plan", type=Path, action="append", default=[])
    b.add_argument("--out", type=Path, required=True)
    b.set_defaults(func=cmd_aggregate)
    return p


if __name__ == "__main__":
    ns = parser().parse_args()
    if getattr(ns, "checkpoint_every", 1) <= 0:
        raise SystemExit("--checkpoint-every must be positive")
    ns.func(ns)
