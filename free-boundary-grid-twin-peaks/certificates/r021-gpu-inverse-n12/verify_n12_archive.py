#!/usr/bin/env python3
"""Cross-check the archived n=12 inverse-fibre finite diagnostic."""
import hashlib
import json
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 90
HERE = Path(__file__).resolve().parent
RAW = HERE / "n12-p65536-r128-t145-l16-domain9300000000000000.jsonl"
SOURCE = HERE / "inverse_fibre_gpu_r020.swift"
SUMMARY = HERE / "n12-p65536-r128-t145-l16.summary.json"
ROOT = HERE / "n12-inverse-root-summary.json"
TELEMETRY = HERE / "n12-p65536-r128-t145-l16.telemetry.jsonl"

EXPECTED_RAW = "13473833ed000fa5737fce1819eff5355b7f4d38ceef3f88b9ca20f7a38769ef"
EXPECTED_SOURCE = "e56997103ba08ea5507b667499f3d059349c6c4b96c00cf1599484bdf78f2e56"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def near(a, b):
    a, b = Decimal(a), Decimal(b)
    return abs(a - b) <= Decimal("1e-14") * max(abs(a), abs(b), Decimal(1))


def main():
    assert sha(RAW) == EXPECTED_RAW
    assert sha(SOURCE) == EXPECTED_SOURCE
    prod = json.loads(SUMMARY.read_text())
    root = json.loads(ROOT.read_text())
    assert prod["raw_sha256"] == EXPECTED_RAW == root["sha256"]
    assert (prod["run"]["n"], prod["run"]["P"], prod["run"]["R"], prod["run"]["tau"]) == (12, 65536, 128, 145)
    assert (root["n"], root["P"], root["R"]) == (12, 65536, 128)
    assert near(prod["probability_normalization"]["mean_Z"], root["mean_Z_onepeak_probability"])
    assert near(prod["probability_normalization"]["A1_estimate_decimal"], root["mean_Z_times_factorial_A1"])
    for d in ("2", "3"):
        assert near(prod["probability_normalization"][f"R{d}_estimate_A2d{d}_over_A1"], root["ratios"][d]["A2_distance_over_A1"])
    records = [json.loads(x) for x in TELEMETRY.read_text().splitlines()]
    exit_row = next(x for x in records if x.get("kind") == "exit")
    assert exit_row["returncode"] == 0
    out = {
        "status": "passed",
        "source_sha256": EXPECTED_SOURCE,
        "raw_sha256": EXPECTED_RAW,
        "production_root_replay_common_values": "equal to relative Decimal tolerance 1e-14",
        "telemetry_exit": exit_row,
        "scope": "finite numerical diagnostic only; no mixing, bias, or asymptotic certificate",
    }
    (HERE / "n12-archive-verified.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
