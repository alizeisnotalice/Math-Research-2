#!/usr/bin/env python3
"""Audit the archived n=7 bridge and n=8 remote receipts without a binary."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "remote-results"
SOURCE = "221ae28c1e074125b2a2fbc9794668c559d3b95f68e6f66de5740304746a6684"
PRIMES = ["9223372036854775783", "9223372036854775643", "9223372036854775549", "9223372036854775507"]
N7 = ["9145417529259655831", "6406629068328171267", "4040678413888520043", "6908390690355257047"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resource(path: Path) -> dict[str, int | float]:
    lines = path.read_text().splitlines()
    timing = next(line for line in lines if " real " in line)
    out: dict[str, int | float] = {"real_seconds": float(timing.split()[0])}
    for line in lines:
        if "maximum resident set size" in line:
            out["max_rss_bytes"] = int(line.split()[0])
        if line.strip().endswith("swaps"):
            out["swaps"] = int(line.split()[0])
        if "peak memory footprint" in line:
            out["peak_memory_footprint_bytes"] = int(line.split()[0])
    return out


def check(tag: str, n: int, index: int) -> dict:
    result_path = RESULTS / f"{tag}.json"
    exit_path = RESULTS / f"{tag}.exit.json"
    start_path = RESULTS / f"{tag}.started.json"
    resource_path = RESULTS / f"{tag}.resource.txt"
    result = json.loads(result_path.read_text())
    exit_data = json.loads(exit_path.read_text())
    start = json.loads(start_path.read_text())
    usage = resource(resource_path)
    checks = {
        "exit_zero": exit_data["exit_code"] == 0,
        "result_hash": sha(result_path) == exit_data["output_sha256"],
        "resource_hash": sha(resource_path) == exit_data["resource_sha256"],
        "command_matches": start["command"] == exit_data["command"],
        "zero_swaps": usage["swaps"] == 0,
        "n": result["n"] == n,
        "prime_index": result["prime_index"] == index,
        "prime": result["prime"] == PRIMES[index],
        "threads": result["threads"] == 12,
        "residue_in_field": 0 <= int(result["residue"]) < int(PRIMES[index]),
    }
    if n == 7:
        checks["bridge_residue"] = result["residue"] == N7[index]
    else:
        checks.update({
            "independent_sets": result["independent_sets"] == 851302029,
            "max_layer_count": result["max_layer_count"] == 162959262,
            "modeled_payload": result["max_modeled_payload_bytes"] == 5803578448,
            "payload_cap": result["payload_limit_bytes"] == 7 * 2**30,
        })
    assert all(checks.values()), {"tag": tag, "checks": checks}
    return {"tag": tag, "result": result, "resource": usage, "checks": checks}


def main() -> None:
    assert sha(HERE.parent / "completion_layer_mod.cpp") == SOURCE
    bridge = [check(f"n7-p{i}", 7, i) for i in range(4)]
    passes = [check(f"n8-p{i}", 8, i) for i in range(4)]
    boundary = json.loads((RESULTS / "boundary-test.json").read_text())
    plan = json.loads((RESULTS / "plan-n8.json").read_text())
    driver = json.loads((RESULTS / "driver.exit.json").read_text())
    launcher = json.loads((RESULTS / "launcher.exit.json").read_text())
    assert boundary == {"N": 64, "FULL_is_all_ones": True, "top_mask": "0xff", "bit63_nonzero": True, "corner63_closed_contains_55_62_63": True}
    assert plan["compute_payload_bytes"] == 5803578448
    assert plan["rebuild_payload_bytes"] == 4524899568
    assert plan["payload_limit_bytes"] == 7 * 2**30
    assert driver["exit_code"] == 0 and launcher["driver_exit_code"] == 0 and launcher["caffeinate_exit_code"] == 0
    out = {
        "all_portable_receipt_checks_pass": True,
        "scope": "Archived receipt audit only; binary is deliberately not archived and CRT is checked separately.",
        "source_sha256": SOURCE,
        "n7_bridge_residues": [row["result"]["residue"] for row in bridge],
        "n8_residues": [row["result"]["residue"] for row in passes],
        "modeled_payload_bytes": 5803578448,
        "max_actual_rss_bytes": max(row["resource"]["max_rss_bytes"] for row in passes),
        "max_actual_peak_memory_footprint_bytes": max(row["resource"]["peak_memory_footprint_bytes"] for row in passes),
        "all_n8_swaps_zero": all(row["resource"]["swaps"] == 0 for row in passes),
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
