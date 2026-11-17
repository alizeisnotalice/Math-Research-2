#!/usr/bin/env python3
"""Independent bookkeeping and exact-decimal audit for the completed 12-window MPS run.

This script deliberately corrects the metadata wording: the final cyclic window
wraps because 12*4096 > 46011.  It proves union coverage, not disjointness.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def case_sha(case: dict) -> str:
    return hashlib.sha256(json.dumps(case, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_verifier(path: Path):
    spec = importlib.util.spec_from_file_location("fraction_verifier", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--fraction-verifier", type=Path, required=True)
    ap.add_argument("--worker", type=Path, required=True)
    ap.add_argument("--prefix-module", type=Path, required=True)
    ap.add_argument("--prefix-plan", type=Path, required=True)
    ap.add_argument("--base-source", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    r = a.run
    config = json.loads((r / "config.json").read_text())
    release = json.loads((r / "release.json").read_text())
    result = json.loads((r / "result.json").read_text())
    checkpoint = json.loads((r / "checkpoint.json").read_text())
    windows_doc = json.loads((r / "window-candidates.json").read_text())
    cpu = json.loads((r / "cpu-prefix-audit.json").read_text())
    plan = json.loads(a.prefix_plan.read_text())
    codes = plan["survivor_codes"]
    B, steps, restart = 4096, 72000, 6000
    assert (config["m"], config["n"], config["starts"], config["steps"], config["restart_every"], config["seed"]) == (96, 144, B, steps, restart, 2026111707)
    assert config["coverage_epochs_required"] == 12 and config["coverage_schedule_complete"] is True
    assert config["prefix_survivor_count"] == len(codes) == 46011
    assert config["prefix_codes_sha256"] == cpu["survivor_codes_sha256"]
    assert sha(a.worker) == config["script_sha256"] == cpu["worker_sha256"]
    assert sha(a.prefix_module) == config["prefix_module_sha256"] == cpu["prefix_module_sha256"]
    assert sha(a.prefix_plan) == config["prefix_plan_sha256"] == cpu["plan_sha256"]
    assert sha(a.base_source) == config["frozen_source_sha256"] == cpu["frozen_source_sha256"]
    assert result["complete"] is True and result["status"] == "complete" and result["last_completed_step"] == steps
    assert release == {"cleanup_called": True, "process_will_exit": True, "terminal_status": "complete", "timed_out": False, "last_completed_step": steps, "remaining_tensor_bytes_before_exit": release["remaining_tensor_bytes_before_exit"], "pid": release["pid"]}
    assert windows_doc["schema"] == "d15-window-candidates-v1" and windows_doc["starts"] == B and windows_doc["restart_every"] == restart
    windows = windows_doc["windows"]
    assert len(windows) == 12
    restart_dir = r / "restarts" if (r / "restarts").is_dir() else r
    restarts = [restart_dir / f"restart_{s:07d}.json" for s in range(restart, steps, restart)]
    assert all(p.is_file() for p in restarts) and len(restarts) == 11
    verifier = load_verifier(a.fraction_verifier)
    pool = set(codes)
    assigned_union, samples, sample_fingerprints, residuals = set(), [], set(), []
    for epoch, window in enumerate(windows):
        assert window["window_epoch"] == epoch
        assert window["window_offset_mod_count"] == (epoch * B) % len(codes)
        assert len(window["candidates"]) == 16
        assert (epoch + 1) * restart - 100 <= window["closing_step"] <= (epoch + 1) * restart
        cohort = [codes[(epoch * B + i) % len(codes)] for i in range(B)]
        assigned_union.update(cohort)
        if epoch < 11:
            snap = json.loads(restarts[epoch].read_text())
            assert snap["old_prefix_epoch"] == epoch and snap["next_prefix_epoch"] == epoch + 1
            assert snap["old_prefix_offset_mod_count"] == window["window_offset_mod_count"]
            assert snap["cases"] == window["candidates"]
        else:
            assert result["cases"] == window["candidates"]
            assert checkpoint["cases"] == window["candidates"]
        for case in window["candidates"]:
            finger = case_sha(case)
            assert finger not in sample_fingerprints
            sample_fingerprints.add(finger)
            assert case["assigned_prefix_code"] in cohort
            residual, _, complement = verifier.validate_case(case, 96, 144, 15, Fraction("0.001"), pool)
            assert residual > 0
            samples.append(case["assigned_prefix_code"])
            residuals.append((residual, complement))
    assert assigned_union == pool
    # The last cyclic window repeats indices 0..3140; this is expected, not a gap.
    repeated_slots = 12 * B - len(pool)
    assert repeated_slots == 3141
    assert len(samples) == len(set(samples)) == 192
    out = {
        "status": "passed_192_saved_exact_decimal_samples_and_cyclic_pool_coverage",
        "scope": "Exact audit of saved JSON decimals, source bindings and cyclic prefix assignment; float/MPS trajectories remain non-certifying numerical exploration.",
        "coverage": {"prefix_pool": len(pool), "windows": 12, "starts_per_window": B, "union_covers_pool": True, "cyclic_repeated_slots": repeated_slots, "windows_are_disjoint": False},
        "samples": {"saved": len(samples), "unique_json_cases": len(sample_fingerprints), "unique_prefix_codes": len(set(samples)), "minimum_exact_decimal_residual": str(float(min(x[0] for x in residuals))), "maximum_designated_complement_error": str(float(max(x[1] for x in residuals)))},
        "bindings": {"config": sha(r / "config.json"), "result": sha(r / "result.json"), "checkpoint": sha(r / "checkpoint.json"), "release": sha(r / "release.json"), "windows": sha(r / "window-candidates.json"), "cpu_audit": sha(r / "cpu-prefix-audit.json"), "worker": sha(a.worker), "prefix_module": sha(a.prefix_module), "prefix_plan": sha(a.prefix_plan), "base_source": sha(a.base_source), "fraction_verifier": sha(a.fraction_verifier), "script": sha(Path(__file__))},
        "release": release,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": out["status"], "samples": len(samples), "prefixes": len(assigned_union), "repeated_slots": repeated_slots}, sort_keys=True))


if __name__ == "__main__":
    main()
