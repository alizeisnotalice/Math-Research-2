#!/usr/bin/env python3
"""Validate a fixed-tau SIS production JSONL and write a compact summary.

The reported ratio is self-normalized and is a numerical diagnostic only.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    raw = args.results.read_bytes()
    rows = [json.loads(x) for x in raw.splitlines()]
    kinds = Counter(x["kind"] for x in rows)
    metadata = next(x for x in rows if x["kind"] == "metadata")
    summary = next(x for x in rows if x["kind"] == "summary")
    n = summary["n"]
    expected = n * n - 2
    assert kinds["metadata"] == 1 and kinds["summary"] == 1
    assert kinds["stratum"] == expected and kinds["trajectory_check"] == expected
    assert {x["tau"] for x in rows if x["kind"] == "stratum"} == set(range(2, n * n))
    assert all(x["ok"] for x in rows if x["kind"] == "trajectory_check")
    batches = [x for x in rows if x["kind"] == "batch"]
    assert len({(x["tau"], x["index"]) for x in batches}) == len(batches)
    out = {
        "status": "complete and internally validated JSONL",
        "result_sha256": hashlib.sha256(raw).hexdigest(),
        "n": n,
        "batch_rows": len(batches),
        "strata": expected,
        "trajectory_checks": expected,
        "metadata": metadata,
        "summary": summary,
        "interpretation": "floating-point deterministic-RNG numerical evidence; no confidence guarantee or proof",
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
