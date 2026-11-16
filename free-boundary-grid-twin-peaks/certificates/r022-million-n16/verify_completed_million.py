#!/usr/bin/env python3
"""Structural and receipt verification for the completed n=16 SMC batch.

This intentionally does not run Metal or estimate a probability.  It verifies
the archived raw rectangle, deterministic seed rule, death layer, monitor
termination receipt, file identities, and the fresh Decimal-80 summary.
"""
import hashlib
import json
import math
from itertools import product
from pathlib import Path


HERE = Path(__file__).resolve().parent
FILES = {
    "source": ("bitset_smc_r017_checked.swift",
               "65aefc323b20995564df7a65656b4bfffbc1570487f8be132fafa916225c2f67"),
    "raw": ("raw-full-n16-p1048576-r16-d7500.jsonl",
            "49ef2fc3ba1b63cf2eacb07ee267685554e37eeac595ded1b0d350bd0af5024a"),
    "telemetry": ("telemetry-full-n16-p1048576-r16-d7500.jsonl",
                  "9a645e3b84a76463ae5160af17a48abf1dc7beb6a34620078c4f9b2f36fd35c9"),
}
METADATA = {
    "kind": "metadata", "implementation": "boundary-bitset-experimental",
    "n": 16, "particles": 1048576, "replicates": 16,
    "tau_lo": 2, "tau_hi": 256, "domain": "7500000000000000",
    "method": "strat", "moment_scale": "U_over_Z",
}
MASK64 = (1 << 64) - 1


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def main():
    paths = {name: HERE / filename for name, (filename, _) in FILES.items()}
    for name, path in paths.items():
        assert sha256(path) == FILES[name][1], f"{name} SHA-256 mismatch"
    source_bytes = paths["source"].read_bytes()
    start = source_bytes.index(b'let metal = """') + len(b'let metal = """')
    end = source_bytes.index(b'"""' + bytes([10]) + b'struct P', start)
    assert hashlib.sha256(source_bytes[start:end]).hexdigest() == (
        "6c59a560e34e746633fa5e6c5ba886c759786d9da2df4d313af9a32258a32181")

    rows = read_jsonl(paths["raw"])
    assert rows[0] == METADATA, "metadata differs from completed batch contract"
    assert rows[-1].get("kind") == "run_end", "raw run_end must be final"
    assert set(rows[-1]) == {"kind", "wall_seconds"}
    assert float(rows[-1]["wall_seconds"]) > 0
    records = rows[1:-1]
    assert len(records) == 4080
    assert all(row.get("kind") == "replicate" for row in records)

    expected = set(product(range(2, 257), range(16)))
    seen = [(int(row["tau"]), int(row["replicate"])) for row in records]
    assert len(seen) == len(set(seen)) == len(expected)
    assert set(seen) == expected
    domain = int(METADATA["domain"], 16)
    dead = []
    for row in records:
        tau, rep = int(row["tau"]), int(row["replicate"])
        seed = (0xA022000000000000 + domain + 16 * 1_000_000
                + tau * 10_000 + rep) & MASK64
        assert row["seed_hex"] == format(seed, "x"), (tau, rep, row["seed_hex"])
        if row["dead"]:
            dead.append((tau, rep))
            assert row["Z"] == 0.0 and row["U"] == 0.0 and row["logZ"] is None
        else:
            assert math.isfinite(float(row["logZ"]))
            assert math.isfinite(float(row["Z"])) and float(row["Z"]) >= 0.0
            assert math.isfinite(float(row["U"]))
    expected_dead = set(product((255, 256), range(16)))
    assert set(dead) == expected_dead and len(dead) == 32

    telemetry = read_jsonl(paths["telemetry"])
    assert telemetry[0].get("kind") == "metadata"
    assert telemetry[0]["memory_scope"] == "global memory_pressure; target ps RSS"
    samples = [row for row in telemetry if row.get("kind") == "sample"]
    assert samples and {int(row["pid"]) for row in samples} == {6859}
    assert telemetry[-1] == {"kind": "exit", "returncode": 0,
                             "epoch_seconds": 1794823107.057085}
    monitor_summary = json.loads((HERE / "monitor-summary.json").read_text())
    assert monitor_summary == {
        "kind": "monitor_summary_from_complete_telemetry",
        "period_seconds": 1,
        "target_pid": 6859,
        "sample_count": len(samples),
        "target_rss_kib_min": min(x["target_rss_kib"] for x in samples),
        "target_rss_kib_max": max(x["target_rss_kib"] for x in samples),
        "system_memory_free_pct_min": min(x["system_memory_free_pct"] for x in samples),
        "system_memory_free_pct_max": max(x["system_memory_free_pct"] for x in samples),
        "exit_returncode": 0,
        "exit_epoch_seconds": 1794823107.057085,
        "scope": "global memory_pressure; target ps RSS",
        "not_measured": ["GPU utilisation", "exclusive GPU ownership"],
    }

    stored = json.loads((HERE / "full-n16-p1048576-r16-d7500-by-replicate.json").read_text())
    fresh = json.loads((HERE / "summary-recomputed.json").read_text())
    assert {k: v for k, v in stored.items() if k != "input"} == {
        k: v for k, v in fresh.items() if k != "input"}, "summary changed beyond input path"
    result = {
        "status": "passed",
        "raw_records": len(records),
        "tau_replicate_rectangle": [2, 256, 16],
        "dead_records": len(dead),
        "dead_taus": [255, 256],
        "run_end_wall_seconds": rows[-1]["wall_seconds"],
        "monitor_pid": 6859,
        "monitor_exit_returncode": 0,
        "telemetry_samples": len(samples),
        "embedded_metal_sha256": "6c59a560e34e746633fa5e6c5ba886c759786d9da2df4d313af9a32258a32181",
        "identities": {name: digest for name, (_, digest) in FILES.items()},
        "stored_and_fresh_summary_equal_except_input_path": True,
        "full_mean_sum_tau_Z_decimal": stored["full_mean_sum_tau_Z_decimal"],
        "full_Z_weighted_U_over_Z": stored["full_Z_weighted_U_over_Z"],
        "full_repeat_logZ_ESS": stored["full_repeat_logZ_ESS"],
    }
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
