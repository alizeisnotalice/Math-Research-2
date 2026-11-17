#!/usr/bin/env python3
"""Verify shipped N=9 checkpoint, coverage, source, gate, and exit evidence.

It reads 326 small JSON checkpoints but never compiles C++ or rescans the
1,327,360 finite hosts.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

TOTALS = {1: 256, 2: 16384, 3: 262144, 4: 1048576}
SOURCE_SHA = "be18e9f8c593648c4e2121305f400db5199ed2fcdde46c4522ddbae27b15b5c4"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    try:
        obj = json.loads(path.read_text())
    except Exception as exc:
        raise SystemExit(f"invalid JSON {path}: {exc}") from exc
    if not isinstance(obj, dict):
        raise SystemExit(f"not an object: {path}")
    return obj


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit("FAIL: " + message)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--json-out", type=Path)
    ns = ap.parse_args()
    root = ns.root.resolve()
    source = root / "source" / "bipartite_prefix_dp.cpp"
    require(digest(source) == SOURCE_SHA, "certified C++ source SHA-256 differs")
    strict = load(root / "evidence" / "strict-aggregate.public.json")
    require(strict.get("status") == "completed-strict-coverage", "aggregate not strict-complete")
    require(strict.get("source_sha256") == SOURCE_SHA, "aggregate source SHA differs")
    paths = sorted((root / "evidence" / "results").glob("shard-*.json"))
    require(len(paths) == 326, f"expected 326 checkpoint files, found {len(paths)}")
    intervals: dict[int, list[tuple[int, int, str]]] = {k: [] for k in TOTALS}
    by_digest: dict[str, str] = {}
    for p in paths:
        row = load(p)
        d = digest(p)
        by_digest[d] = p.name
        left, start, stop, next_mask, scanned = (row.get(k) for k in ("left_size", "range_start", "range_stop_exclusive", "next_mask", "graphs_scanned"))
        require(row.get("status") == "completed" and row.get("mode") == "shard-r27" and row.get("N") == 9, f"bad status/mode/N in {p.name}")
        require(row.get("source_sha256") == SOURCE_SHA, f"source mismatch in {p.name}")
        require(type(left) is int and left in TOTALS and all(type(v) is int for v in (start, stop, next_mask, scanned)), f"integer fields malformed in {p.name}")
        require(0 <= start <= next_mask == stop <= TOTALS[left], f"interval bounds invalid in {p.name}")
        require(scanned == stop - start and row.get("max_abs_R27") == 0 and "counterexample" not in row, f"not a certified zero interval in {p.name}")
        intervals[left].append((start, stop, d))
    expected_coverage: dict[str, list[dict]] = {}
    for left, total in TOTALS.items():
        rows = sorted(intervals[left])
        cursor = 0
        rendered: list[dict] = []
        for start, stop, d in rows:
            require(start == cursor, f"left {left} coverage gap/overlap at {cursor}")
            rendered.append({"start": start, "stop_exclusive": stop, "origin": "completed:" + d})
            cursor = stop
        require(cursor == total, f"left {left} has only [0,{cursor}) of [0,{total})")
        expected_coverage[str(left)] = rendered
    require(strict.get("coverage") == expected_coverage, "shipped strict aggregate does not match checkpoint coverage")
    require(strict.get("checkpoint_count") == 326, "aggregate checkpoint count differs")

    exits = load(root / "evidence" / "execution-exits.public.json")
    require(exits.get("source_sha256") == SOURCE_SHA and exits.get("halted_left_sizes") == [], "exit ledger source/halt condition fails")
    all_exits = [exits.get("benchmark")] + exits.get("production", [])
    require(len(all_exits) == 326, f"expected 326 exit records, found {len(all_exits)}")
    ledger_digests: set[str] = set()
    for row in all_exits:
        require(isinstance(row, dict) and row.get("exit_code") == 0 and row.get("checkpoint_status") == "completed", "nonzero/noncompleted execution record")
        name = Path(row.get("checkpoint", "")).name
        p = root / "evidence" / "results" / name
        require(p.is_file() and digest(p) == row.get("checkpoint_sha256"), f"exit record checkpoint mismatch: {name}")
        ledger_digests.add(row["checkpoint_sha256"])
    require(ledger_digests == set(by_digest), "exit ledger is not a one-to-one match to results")

    receipts = load(root / "evidence" / "run-receipts.public.json")
    require(receipts.get("source_sha256") == SOURCE_SHA, "receipt source differs")
    require(receipts.get("compile", {}).get("exit_code") == 0 and receipts.get("aggregate", {}).get("exit_code") == 0, "compile/aggregate exit nonzero")
    require(receipts.get("completed", {}).get("status") == "completed-strict-coverage", "completion status differs")
    gate_result_digests = {"gate-n7-r45": "82f2bd8b654035224d51467cd10f857ff1b459b641a5f6ed89ff775d66fe971e", "gate-n9-direct": "a557a75a580b416e9d3a1c4301216b2e22e3d6aba14763ff7641b73696da4507"}
    for mode, expected in gate_result_digests.items():
        p = root / "evidence" / "gates" / f"{mode}.json"
        row = load(p)
        require(digest(p) == expected and row.get("status") == "passed" and row.get("source_sha256") == SOURCE_SHA, f"gate result invalid: {mode}")
    gates = {g.get("mode"): g for g in receipts.get("gates", [])}
    require(set(gates) == set(gate_result_digests), "gate receipt set differs")
    require(all(g.get("exit_code") == 0 for g in gates.values()), "gate exit nonzero")
    direct = load(root / "evidence" / "gates" / "gate-n9-direct.json")
    samples = direct.get("samples", [])
    require(len(samples) == 3 and all(x.get("dp_equals_direct") is True and x.get("R27") == 0 for x in samples), "three direct-9! samples do not all agree")
    result = {"status": "passed", "certified_source_sha256": SOURCE_SHA, "checkpoint_count": len(paths), "covered_edge_subsets": sum(TOTALS.values()), "coverage_counts": {str(k): len(intervals[k]) for k in TOTALS}, "execution_records": len(all_exits), "n9_direct_samples": [x["name"] for x in samples]}
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if ns.json_out:
        ns.json_out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
