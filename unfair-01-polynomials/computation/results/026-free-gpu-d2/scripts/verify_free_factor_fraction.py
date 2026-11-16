#!/usr/bin/env python3
"""Independent exact-rational audit of saved free-coefficient GPU samples.

No worker/search module is imported.  JSON decimal coefficients are treated as
their exact displayed rational values; this validates saved records, not the
unrecorded float32 tensors or an exclusion/counterexample claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(x) -> Fraction:
    if not isinstance(x, (int, float)):
        raise ValueError(f"coefficient is not JSON number: {x!r}")
    return Fraction(str(x))


def convolution(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    c = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def audit_case(case: dict, m: int, n: int, d: int, margin: Fraction) -> dict:
    a, b = [q(x) for x in case["A"]], [q(x) for x in case["B"]]
    if len(a) != m + 1 or len(b) != n + 1:
        raise ValueError("factor degree/length mismatch")
    if a[0] != 1 or a[-1] != 1 or b[0] != 1 or b[-1] != 1:
        raise ValueError("endpoint model mismatch")
    if any(x < 0 or x > 1 for x in a + b):
        raise ValueError("coefficient leaves [0,1]")
    if any(x != 0 for x in a[1:d] + b[1:d]):
        raise ValueError("forced low-prefix zero mismatch")
    if not (margin <= a[d] <= 1 - margin and margin <= b[d] <= 1 - margin):
        raise ValueError("designated coefficient violates margin")
    complement_delta = abs(a[d] + b[d] - 1)
    if complement_delta > Fraction(1, 10**7):
        raise ValueError("designated complementary pair differs from one by more than float32 tolerance")
    c = convolution(a, b)
    target = [Fraction(int(x >= Fraction(1, 2))) for x in c]
    residual = max(abs(x - y) for x, y in zip(c, target))
    stored = q(case["max_residual_float64"])
    difference = abs(residual - stored)
    # JSON has a rounded binary64 diagnostic, whereas the calculation above is
    # over the exact displayed decimal values.  This is a consistency gate, not
    # an assertion that the two arithmetic models are identical.
    if difference > Fraction(1, 10**12):
        raise ValueError(f"stored max residual differs by {float(difference):.3g}")
    binary = all(x in (0, 1) for x in a + b)
    fractional = any(0 < x < 1 for x in a + b)
    return {
        "max_residual_fraction": f"{residual.numerator}/{residual.denominator}",
        "max_residual_decimal": str(float(residual)),
        "stored_residual_delta": str(float(difference)),
        "exact_binary_factors": binary,
        "has_fractional_coefficient": fractional,
        "designated_A_d": str(float(a[d])),
        "designated_B_d": str(float(b[d])),
        "designated_complement_delta": str(float(complement_delta)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--record", type=Path, required=True)
    args = ap.parse_args()
    run = args.run_dir
    config_path, checkpoint_path = run / "config.json", run / "checkpoint.json"
    config, checkpoint = json.loads(config_path.read_text()), json.loads(checkpoint_path.read_text())
    m, n, d = (int(config[k]) for k in ("m", "n", "fractional_index"))
    margin = q(config["margin"])
    restart_dir = run / "restarts" if (run / "restarts").is_dir() else run
    restart_paths = sorted(restart_dir.glob("restart_*.json"))
    expected = [restart_dir / f"restart_{i:07d}.json" for i in range(int(config["restart_every"]), int(config["steps"]), int(config["restart_every"]))]
    if restart_paths != expected:
        raise SystemExit("restart ledger does not have exactly the configured pre-final rotation snapshots")
    docs = [(p.name, json.loads(p.read_text()), "complete_rotation_snapshot") for p in restart_paths]
    docs.append(("checkpoint.json", checkpoint, "complete_final_cohort"))
    summaries, seen = [], set()
    rotation_rows = final_rows = 0
    max_saved_delta = Fraction(0)
    max_complement_delta = Fraction(0)
    exact_binary_count = fractional_count = 0
    for name, doc, kind in docs:
        cases = doc.get("cases")
        hist = doc.get("history")
        if not isinstance(cases, list) or len(cases) != 16 or not isinstance(hist, list) or not hist:
            raise SystemExit(f"{name}: expected exactly 16 cases and nonempty history")
        if kind == "complete_rotation_snapshot" and len(hist) != 1:
            raise SystemExit(f"{name}: restart snapshot must contain one local history row")
        row_results = [audit_case(case, m, n, d, margin) for case in cases]
        fingerprints = [hashlib.sha256(json.dumps(case, sort_keys=True, separators=(",", ":")).encode()).hexdigest() for case in cases]
        if len(set(fingerprints)) != 16:
            raise SystemExit(f"{name}: duplicate saved case")
        seen.update(fingerprints)
        max_saved_delta = max(max_saved_delta, *(Fraction(x["stored_residual_delta"]) for x in row_results))
        max_complement_delta = max(max_complement_delta, *(Fraction(x["designated_complement_delta"]) for x in row_results))
        exact_binary_count += sum(x["exact_binary_factors"] for x in row_results)
        fractional_count += sum(x["has_fractional_coefficient"] for x in row_results)
        if kind == "complete_rotation_snapshot": rotation_rows += len(cases)
        else: final_rows += len(cases)
        summaries.append({"file": name, "kind": kind, "history_step": hist[-1].get("step"),
                          "cases": len(cases), "minimum_fraction_residual": min(x["max_residual_decimal"] for x in row_results),
                          "maximum_stored_residual_delta": str(float(max(Fraction(x["stored_residual_delta"]) for x in row_results)))})
    release_path = run / "release.json"
    release = json.loads(release_path.read_text())
    if release.get("cleanup_called") is not True or release.get("process_will_exit") is not True:
        raise SystemExit("release record does not attest cleanup/exit")
    result = {
        "status": "saved_samples_fraction_reconvolution_passed",
        "claim": "Saved JSON sample audit only; numeric free-coefficient search is non-strict and does not prove exclusion or provide an exact counterexample.",
        "run_hashes": {"config": digest(config_path), "checkpoint": digest(checkpoint_path), "release": digest(release_path)},
        "model": {"m": m, "n": n, "fractional_index": d, "margin": str(float(margin)), "starts_per_cohort": int(config["starts"]), "restart_every": int(config["restart_every"])},
        "counts": {"restart_files": len(restart_paths), "rotation_snapshot_saved_rows": rotation_rows,
                   "complete_final_cohort_saved_rows": final_rows, "partial_saved_rows": 0, "total_saved_rows": rotation_rows + final_rows,
                   "unique_saved_json_cases": len(seen), "exact_binary_factor_pairs": exact_binary_count,
                   "pairs_with_a_fractional_coefficient": fractional_count},
        "max_stored_residual_delta": str(float(max_saved_delta)),
        "max_designated_complement_delta": str(float(max_complement_delta)),
        "release_before_exit": {"cleanup_called": release["cleanup_called"], "remaining_tensor_bytes_before_exit": release["remaining_tensor_bytes_before_exit"],
                                "process_will_exit": release["process_will_exit"], "audit_note": "GPU process was independently absent at audit time; this nonzero value was recorded before exit."},
        "files": summaries,
    }
    args.record.parent.mkdir(parents=True, exist_ok=True)
    args.record.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
