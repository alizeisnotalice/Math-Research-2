#!/usr/bin/env python3
"""MP80 audit for all retained DS rows, including missing CPU copies.

This deliberately supplements, rather than relaxes, verify_coprime_result.py:
the latter still rejects a checkpoint with missing CPU copies.  We
independently decode the same explicit option mapping and reconstruct every
row at MP precision, comparing stored copies when they exist.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import mpmath as mp


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_verifier(path: Path):
    spec = importlib.util.spec_from_file_location("coprime_result_verifier", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--options", type=Path, required=True)
    ap.add_argument("--roots-input", type=Path, required=True)
    ap.add_argument("--engine-source", type=Path, required=True)
    ap.add_argument("--verifier-source", type=Path, required=True)
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--mp-dps", type=int, default=80)
    ap.add_argument("--tolerance", type=str, default="1e-7")
    args = ap.parse_args()
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data not mounted")
    if args.mp_dps < 50:
        raise SystemExit("MP precision below 50 refused")
    mp.mp.dps = args.mp_dps
    tol = mp.mpf(args.tolerance)
    mod = load_verifier(args.verifier_source)
    checkpoint = json.loads(args.checkpoint.read_text())
    options = json.loads(args.options.read_text())
    raw = json.loads(args.roots_input.read_text())["candidates"][int(checkpoint["candidate_index"])]
    # Written without a broad 'success' gate: every identity is checked here.
    if checkpoint.get("format") != "coprime-substitution-double-single-mps-v1":
        raise SystemExit("unexpected checkpoint format")
    if checkpoint.get("source_sha256") != sha(args.engine_source):
        raise SystemExit("engine hash mismatch")
    if checkpoint.get("options_sha256") != sha(args.options):
        raise SystemExit("options hash mismatch")
    if checkpoint.get("roots_input_sha256") != sha(args.roots_input):
        raise SystemExit("roots-input hash mismatch")
    if checkpoint.get("candidate") != raw["name"]:
        raise SystemExit("candidate mismatch")
    a, b = int(raw["a"]), int(raw["b"])
    layouts = mod.layouts(options, a, b)
    eta = mp.mpf(str(checkpoint["parameters"]["eta"]))
    target = [mp.mpf(x) for x in raw["C_coefficients_ascending"]]
    records = []
    max_residual = mp.mpf("0")
    for entry in checkpoint.get("E_progress", []):
        E = int(entry["E"])
        radices, total = layouts[E]
        coverage = entry["coverage_completed_half_open"]
        if coverage[0] != 0 or coverage[1] > total or entry["total_codes"] != total or entry["radices_block0_first"] != radices:
            raise SystemExit(f"bad ledger at E={E}")
        for row in entry.get("top64", []):
            code = int(row["code"])
            if not 0 <= code < coverage[1]:
                raise SystemExit(f"row outside prefix E={E}, code={code}")
            digits = mod.decode(code, radices)
            if [int(x) for x in row["mixed_radix_digits_block0_first"]] != digits:
                raise SystemExit(f"decode mismatch E={E}, code={code}")
            selected, mapping = mod.expected_mapping(options, E, digits)
            if mapping != row["mask_mapping"]:
                raise SystemExit(f"mapping mismatch E={E}, code={code}")
            bitmask = sum(1 << index for index in selected)
            if row["A_root_group_bitmask_decimal"] != str(bitmask):
                raise SystemExit(f"mask mismatch E={E}, code={code}")
            left, right, endpoints = mod.reconstruct(raw["groups"], selected)
            saved = row.get("cpu_float64_independent_root_reconstruction")
            saved_difference = None
            saved_fractionality_match = None
            if saved:
                saved_difference = max(
                    mod.max_difference(left, saved["A_coefficients_ascending"]),
                    mod.max_difference(right, saved["B_coefficients_ascending"]),
                )
                if saved_difference > tol:
                    raise SystemExit(f"MP saved-polynomial mismatch {saved_difference} at E={E}, code={code}")
            product = mod.poly_mul(left, right)
            residual = max(abs(x - y) for x, y in zip(product, target))
            if residual > tol:
                raise SystemExit(f"MP product residual {residual} exceeds tolerance at E={E}, code={code}")
            metrics = mod.mp_metrics(left, right, eta)
            if saved and bool(saved["has_fractionality_eta_1e-5"]) != bool(metrics["fractional"]):
                raise SystemExit(f"saved fractionality mismatch E={E}, code={code}")
            max_residual = max(max_residual, residual)
            records.append({
                "E": E, "code": code, "coverage_prefix_end": coverage[1],
                "fractional_eta": bool(metrics["fractional"]),
                "minimum": str(metrics["minimum"]), "maximum": str(metrics["maximum"]),
                "score": str(metrics["score"]), "convolution_residual": str(residual),
                "stored_cpu_copy_present": bool(saved),
                "stored_cpu_max_difference": None if saved_difference is None else str(saved_difference),
                "raw_endpoint_deviation": [str(x) for x in endpoints],
            })
    if not records:
        raise SystemExit("no missing stored-CPU rows found")
    out = {
        "status": "all_retained_rows_reconstructed_mp",
        "claim_boundary": "Supplements missing stored-CPU fields with MP reconstruction; it neither changes the frozen checkpoint nor turns a partial prefix into full coverage.",
        "checkpoint": str(args.checkpoint), "checkpoint_sha256": sha(args.checkpoint),
        "engine_source_sha256": sha(args.engine_source), "options_sha256": sha(args.options),
        "roots_input_sha256": sha(args.roots_input), "verifier_source_sha256": sha(args.verifier_source),
        "script_sha256": sha(Path(__file__)), "a": a, "b": b, "mp_dps": args.mp_dps,
        "tolerance": str(tol), "rows": records, "row_count": len(records),
        "fractional_rows": sum(x["fractional_eta"] for x in records),
        "stored_cpu_rows": sum(x["stored_cpu_copy_present"] for x in records),
        "missing_stored_cpu_rows": sum(not x["stored_cpu_copy_present"] for x in records),
        "max_convolution_residual": str(max_residual), "gpu_used": False,
    }
    args.record.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.record.with_suffix(args.record.suffix + ".tmp")
    tmp.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.record)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
