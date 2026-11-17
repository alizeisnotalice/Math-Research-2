#!/usr/bin/env python3
"""Decompress the public n=32 raw, invoke frozen f1593 auditor, compare fields."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path

FROZEN_AUDITOR_SHA256 = "f1593cc93db21d778289966c80539606b079d702fa3fd5c197394c27d4737389"
FIELDS = ("n", "P", "R", "L", "ESS", "maxshare", "m_lift_ESS", "m_lift_maxshare",
          "scaled_mean_m", "scaled_mean_mD", "scaled_mean_mF", "scaled_mean_mG",
          "A2_over_B_mRB", "Q_mRB", "Q_L_legacy", "log_B_over_Nfactorial",
          "log_mA2_over_Nfactorial")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a: object, b: object) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(float(a), float(b), rel_tol=2e-13, abs_tol=2e-13)
    return a == b


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-gz", type=Path, required=True)
    parser.add_argument("--auditor", type=Path, required=True)
    parser.add_argument("--derived-source", type=Path, required=True,
                        help="public historical Swift source whose fingerprint the frozen auditor requires")
    parser.add_argument("--recorded", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert sha256(args.auditor) == FROZEN_AUDITOR_SHA256, "frozen auditor SHA-256 mismatch"
    args.workdir.mkdir(parents=True, exist_ok=True)
    raw = args.workdir / "n32-p262144-r32-l16.jsonl"
    with gzip.open(args.raw_gz, "rb") as src, raw.open("wb") as dst:
        dst.write(src.read())
    receipt = json.loads(args.receipt.read_text())
    assert sha256(raw) == receipt["raw_sha256"], "decompressed raw SHA-256 mismatch"
    spec = importlib.util.spec_from_file_location("frozen_auditor", args.auditor)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    replay = module.audit(raw, args.derived_source)
    recorded = json.loads(args.recorded.read_text())
    mismatches = {key: [replay.get(key), recorded.get(key)] for key in FIELDS
                  if not close(replay.get(key), recorded.get(key))}
    assert not mismatches, mismatches
    result = {
        "status": "passed",
        "scope": "raw integrity and frozen finite implementation-invariant audit only; no cone check, probability, mixing, or SMC-accuracy certification",
        "auditor_sha256": FROZEN_AUDITOR_SHA256,
        "raw_sha256": sha256(raw),
        "fields_compared": list(FIELDS),
        "field_comparison": "passed",
        "replay": replay,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "raw_sha256": result["raw_sha256"],
                      "field_comparison": result["field_comparison"]}, indent=2))


if __name__ == "__main__":
    main()
