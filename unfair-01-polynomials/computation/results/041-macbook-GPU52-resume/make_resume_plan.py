#!/usr/bin/env python3
"""Build a restart input without silently crediting an interrupted stream.

This program is intentionally local-only.  It checks the received old GPU52
batch stream before assigning the first unreceived mask of each task.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

EXPECTED_STARTS = {0: 67108864, 1: 67108864, 2: 23166976}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_old_prefix(stdout: Path, tasks: list[dict]) -> dict[int, int]:
    by_id = {int(t["id"]): int(t["masks"]) for t in tasks}
    cursor = {i: 0 for i in by_id}
    for line_no, line in enumerate(stdout.read_text().splitlines(), 1):
        # Wolfram's GPU compatibility progress is deliberately not evidence.
        # The prior independent audit established that only JSON event rows
        # participate in the received-prefix ledger.
        if not line.startswith("{"):
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError as exc:
            raise SystemExit(f"old stdout line {line_no}: invalid JSON") from exc
        if row.get("event") != "batch":
            continue
        task = int(row["task"])
        if task not in by_id:
            raise SystemExit(f"old stdout line {line_no}: unknown task {task}")
        first, last, count = int(row["first"]), int(row["last"]), int(row["count"])
        if first != cursor[task] or last - first + 1 != count:
            raise SystemExit(
                f"old stdout line {line_no}: non-contiguous task {task}: "
                f"expected {cursor[task]}, got {first}..{last}"
            )
        cursor[task] = last + 1
        if cursor[task] > by_id[task]:
            raise SystemExit(f"old stdout line {line_no}: task {task} exceeds masks")
    if {i: cursor[i] for i in EXPECTED_STARTS} != EXPECTED_STARTS:
        raise SystemExit(f"old received prefix differs from audited prefix: {cursor}")
    if any(cursor[i] != 0 for i in by_id if i not in EXPECTED_STARTS):
        raise SystemExit("old stream unexpectedly contains a later task")
    return cursor


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--old-stdout", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--receipt", type=Path, required=True)
    args = ap.parse_args()
    raw = json.loads(args.original.read_text())
    tasks = raw.get("tasks")
    if not isinstance(tasks, list) or [t.get("id") for t in tasks] != list(range(24)):
        raise SystemExit("expected the original 24 tasks in ids 0..23")
    cursor = read_old_prefix(args.old_stdout, tasks)
    resumed = dict(raw)
    resumed["tasks"] = [dict(t, resume_start_mask=cursor[int(t["id"])]) for t in tasks]
    # This source binding is metadata, not a substitute for independently
    # checking the old stream (which this program just did).
    resumed["resume_provenance"] = {
        "original_input_sha256": sha256(args.original),
        "old_stdout_sha256": sha256(args.old_stdout),
        "semantics": "Only batch records in the old stdout establish received prefix; skipped full tasks are not recomputed.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(resumed, indent=2, sort_keys=True) + "\n")
    receipt = {
        "event": "resume_plan_built",
        "original": str(args.original),
        "original_sha256": sha256(args.original),
        "old_stdout": str(args.old_stdout),
        "old_stdout_sha256": sha256(args.old_stdout),
        "resume_input": str(args.out),
        "resume_input_sha256": sha256(args.out),
        "starts": {str(i): cursor[i] for i in range(24)},
        "covered_old_masks": sum(cursor.values()),
        "remaining_masks": sum(int(t["masks"]) - cursor[int(t["id"])] for t in tasks),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
