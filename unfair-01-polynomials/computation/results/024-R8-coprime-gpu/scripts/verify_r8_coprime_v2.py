#!/usr/bin/env python3
"""Independent, CPU-only audit for saved R8 four-block v2 scan records.

This file intentionally does not import the GPU search program.  It checks the
v2 JSON ledger, reconstructs every retained top-64 sample from the frozen root
input using mpmath, and labels an incomplete prefix as such.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import mpmath as mp


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def franel(order: int, blocks: int) -> int:
    return sum(math.comb(order, k) ** blocks for k in range(order + 1))


def expected_layout_total(a: int, b: int, blocks: int) -> int:
    # Two terminal allocations have 2^(blocks-1) canonical layouts, whereas
    # their contribution to F_blocks(a)F_blocks(b)/2 is two.  Hence +2^(n-1)-2.
    return int(franel(a, blocks) * franel(b, blocks) // 2 + (1 << (blocks - 1)) - 2)


def conv(p: list[mp.mpf], q: list[mp.mpf]) -> list[mp.mpf]:
    out = [mp.mpf("0")] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def quadratic_from_lift_root(group: dict) -> list[mp.mpf]:
    """Reconstruct the monic real quadratic from the stored high-precision root."""
    root = group.get("lift_root")
    if not isinstance(root, list) or len(root) != 2:
        raise ValueError(f"group {group.get('index')} has no complex lift_root pair")
    re, im = mp.mpf(str(root[0])), mp.mpf(str(root[1]))
    return [re * re + im * im, -2 * re, mp.mpf("1")]


def selected_polynomial(groups: dict[int, dict], indices: list[int]) -> list[mp.mpf]:
    out = [mp.mpf("1")]
    for index in indices:
        if index not in groups:
            raise ValueError(f"unknown group index {index}")
        out = conv(out, quadratic_from_lift_root(groups[index]))
    return out


def reconstruct(groups: dict[int, dict], a_indices: list[int], all_indices: list[int]):
    a_poly = selected_polynomial(groups, a_indices)
    a_set = set(a_indices)
    b_poly = selected_polynomial(groups, [i for i in all_indices if i not in a_set])
    a0, b0 = a_poly[0], b_poly[0]
    a_poly = [x / a0 for x in a_poly]
    b_poly = [x / b0 for x in b_poly]
    endpoint_deviation = max(abs(a_poly[-1] - 1), abs(b_poly[-1] - 1))
    return a_poly, b_poly, endpoint_deviation


def decode(code: int, radices: list[int]) -> list[int]:
    digits = []
    for radix in radices:
        digits.append(code % radix)
        code //= radix
    if code:
        raise ValueError("code exceeds mixed-radix range")
    return digits


def mask(indices: list[int]) -> int:
    return sum(1 << i for i in indices)


def metric(a_poly: list[mp.mpf], b_poly: list[mp.mpf], eta: mp.mpf):
    interior = a_poly[1:-1] + b_poly[1:-1]
    actual = a_poly + b_poly
    if not actual:
        return mp.mpf("0"), mp.mpf("0"), mp.mpf("0"), False
    minimum, maximum = min(actual), max(actual)
    score = sum(max(-x, mp.mpf("0")) + max(x - 1, mp.mpf("0")) for x in interior)
    fractional = any(eta < x < 1 - eta for x in interior)
    return minimum, maximum, score, fractional


def cpu_copy_compare(saved: dict, a_poly: list[mp.mpf], b_poly: list[mp.mpf], residual: mp.mpf,
                     minimum: mp.mpf, maximum: mp.mpf, fractional: bool, tol: mp.mpf) -> list[str]:
    errs: list[str] = []
    pairs = (("A_coefficients_ascending", a_poly), ("B_coefficients_ascending", b_poly))
    for key, now in pairs:
        old = saved.get(key)
        if not isinstance(old, list) or len(old) != len(now):
            errs.append(f"stored CPU {key} malformed")
            continue
        delta = max(abs(mp.mpf(str(x)) - y) for x, y in zip(old, now))
        if delta > tol:
            errs.append(f"stored CPU {key} delta {mp.nstr(delta, 8)} > {mp.nstr(tol, 4)}")
    scalar_pairs = (
        ("convolution_max_abs_residual", residual),
        ("minimum_actual_coefficient", minimum),
        ("maximum_actual_coefficient", maximum),
    )
    for key, now in scalar_pairs:
        if key not in saved:
            errs.append(f"stored CPU missing {key}")
        elif abs(mp.mpf(str(saved[key])) - now) > tol:
            errs.append(f"stored CPU {key} mismatch")
    if saved.get("has_fractionality_eta_1e-5") is not fractional:
        errs.append("stored CPU fractionality flag mismatch")
    return errs


def check_layouts(options: dict, a: int, b: int, blocks: int, terminal_E: int, errors: list[str]):
    block_data = options.get("blocks")
    if not isinstance(block_data, list) or len(block_data) != blocks:
        fail(errors, f"options blocks has wrong length: {len(block_data) if isinstance(block_data, list) else 'not-list'}")
        return {}, {}
    layouts, all_group_indices = {}, []
    e_keys = None
    for ordinal, block in enumerate(block_data):
        block_id = block.get("block")
        if block_id != ordinal:
            fail(errors, f"noncanonical/missing block index {block_id}")
        all_ids = block.get("all_group_indices", [])
        a_ids, b_ids = block.get("a_lift_groups", []), block.get("b_lift_groups", [])
        if set(all_ids) != set(a_ids) | set(b_ids) or set(a_ids) & set(b_ids):
            fail(errors, f"block {block_id} family partition fails")
        all_group_indices.extend(all_ids)
        table = block.get("options", {})
        if e_keys is None:
            e_keys = set(table)
        elif set(table) != e_keys:
            fail(errors, f"block {block_id} E-key mismatch")
        for e_text, choices in table.items():
            E = int(e_text)
            if not choices:
                fail(errors, f"E={E}, block={block_id}: no choices")
                continue
            terminal = E == terminal_E
            if terminal:
                if {x.get("mode") for x in choices} != {"all_a", "all_b"} or len(choices) != 2:
                    fail(errors, f"E={E}, block={block_id}: malformed terminal choices")
            else:
                for choice in choices:
                    if choice.get("mode") != "mixed" or int(choice.get("u", -1)) * b + int(choice.get("v", -1)) * a != E:
                        fail(errors, f"E={E}, block={block_id}: invalid balance option")
            layouts.setdefault(E, []).append(choices)
    if sorted(all_group_indices) != list(range(len(all_group_indices))):
        fail(errors, "global group indices are not exactly contiguous")
    return layouts, {i for i in all_group_indices}


def audit(args) -> tuple[dict, int]:
    paths = {name: Path(getattr(args, name)) for name in ("checkpoint", "options", "roots_input", "engine_source")}
    errors: list[str] = []
    for name, path in paths.items():
        if not path.is_file():
            fail(errors, f"missing {name}: {path}")
    if errors:
        return {"result": "FAIL", "errors": errors}, 2
    checkpoint = json.loads(paths["checkpoint"].read_text())
    options = json.loads(paths["options"].read_text())
    root_data = json.loads(paths["roots_input"].read_text())
    hashes = {name: sha256_file(path) for name, path in paths.items()}
    a, b, blocks = int(options["a"]), int(options["b"]), int(options["n_blocks"])
    terminal_E = int(options["balance_terminal_E"])
    if blocks != 4:
        fail(errors, f"this R8 audit requires n_blocks=4, found {blocks}")
    if checkpoint.get("source_sha256") != hashes["engine_source"]:
        fail(errors, "checkpoint engine source SHA does not agree with frozen source")
    if checkpoint.get("options_sha256") != hashes["options"]:
        fail(errors, "checkpoint options SHA mismatch")
    if checkpoint.get("roots_input_sha256") != hashes["roots_input"]:
        fail(errors, "checkpoint roots-input SHA mismatch")
    layouts, expected_ids = check_layouts(options, a, b, blocks, terminal_E, errors)
    candidate_name = options.get("candidate")
    candidates = [x for x in root_data.get("candidates", []) if x.get("name") == candidate_name]
    if len(candidates) != 1:
        fail(errors, f"expected exactly one root candidate {candidate_name!r}, found {len(candidates)}")
        candidate = {}
    else:
        candidate = candidates[0]
    groups = {int(x["index"]): x for x in candidate.get("groups", [])}
    if set(groups) != expected_ids:
        fail(errors, "root-input group index set does not agree with options")
    expected_group_count = blocks * (a + b)
    if len(groups) != expected_group_count:
        fail(errors, f"expected {expected_group_count} groups from n_blocks*(a+b), found {len(groups)}")
    target = [mp.mpf(x) for x in candidate.get("C_coefficients_ascending", [])]
    if not target or any(x not in (0, 1) for x in target):
        fail(errors, "target C is not a nonempty {0,1} coefficient list")
    mp.mp.dps = args.mp_dps
    tol, eta = mp.mpf(str(args.tolerance)), mp.mpf("1e-5")
    max_lift_root_float64_delta = mp.mpf("0")
    for group in groups.values():
        try:
            exact_quad = quadratic_from_lift_root(group)
            replica = group.get("coefficients_ascending")
            if not isinstance(replica, list) or len(replica) != 3:
                fail(errors, f"group {group.get('index')}: malformed float64 quadratic replica")
                continue
            delta = max(abs(x - mp.mpf(str(y))) for x, y in zip(exact_quad, replica))
            max_lift_root_float64_delta = max(max_lift_root_float64_delta, delta)
            if delta > mp.mpf("1e-12"):
                fail(errors, f"group {group.get('index')}: lift-root/float64 quadratic delta {mp.nstr(delta, 8)} > 1e-12")
        except Exception as exc:
            fail(errors, f"group {group.get('index')}: lift-root reconstruction exception {type(exc).__name__}: {exc}")
    rows = 0
    missing_cpu_rows = 0
    cpu_copy_rows = 0
    max_residual = mp.mpf("0")
    max_endpoint_deviation = mp.mpf("0")
    closed_rows = 0
    per_E = []
    E_progress = checkpoint.get("E_progress")
    if not isinstance(E_progress, list):
        fail(errors, "checkpoint E_progress is not a list")
        E_progress = []
    for progress in E_progress:
        E = int(progress.get("E", -1))
        choices_by_block = layouts.get(E)
        if choices_by_block is None:
            fail(errors, f"checkpoint E={E} missing from options")
            continue
        radices = [len(x) for x in choices_by_block]
        terminal = E == terminal_E
        expected_total = 1
        expected_radices = radices[:]
        if terminal:
            expected_radices[0] = 1
        for radix in expected_radices:
            expected_total *= radix
        if progress.get("radices_block0_first") != expected_radices:
            fail(errors, f"E={E}: stored radices disagree with option table")
        if int(progress.get("total_codes", -1)) != expected_total:
            fail(errors, f"E={E}: total_codes mismatch")
        coverage = progress.get("coverage_completed_half_open")
        if not (isinstance(coverage, list) and len(coverage) == 2 and coverage[0] == 0 and 0 <= coverage[1] <= expected_total):
            fail(errors, f"E={E}: invalid coverage prefix")
            coverage = [0, 0]
        counts = progress.get("counts", {})
        if int(counts.get("codes_seen", -1)) != coverage[1]:
            fail(errors, f"E={E}: codes_seen does not equal coverage end")
        top = progress.get("top64", [])
        if not isinstance(top, list) or len(top) > 64:
            fail(errors, f"E={E}: malformed top64")
            top = []
        per_E.append({"E": E, "status": progress.get("status"), "coverage": coverage,
                      "total_codes": expected_total, "saved_rows": len(top),
                      "stored_cpu_copies": sum("cpu_float64_independent_root_reconstruction" in x for x in top)})
        for item in top:
            rows += 1
            code = int(item.get("code", -1))
            if not 0 <= code < coverage[1]:
                fail(errors, f"E={E} code={code}: not inside saved coverage prefix")
                continue
            try:
                digits = decode(code, expected_radices)
            except ValueError as exc:
                fail(errors, f"E={E} code={code}: {exc}")
                continue
            selected = []
            expected_map = []
            for block, digit in enumerate(digits):
                choice = choices_by_block[block][digit]
                chosen = list(choice.get("selected_group_indices", []))
                selected.extend(chosen)
                expected_map.append({"block": block, "option_index": digit, "mode": choice.get("mode"),
                                     "A_group_indices": chosen})
            if item.get("mixed_radix_digits_block0_first") != digits:
                fail(errors, f"E={E} code={code}: digit mismatch")
            if item.get("mask_mapping") != expected_map:
                fail(errors, f"E={E} code={code}: option mapping mismatch")
            actual_mask = mask(selected)
            if int(item.get("A_root_group_bitmask_decimal", -1)) != actual_mask:
                fail(errors, f"E={E} code={code}: decimal bitmask mismatch")
            if item.get("A_root_group_bitmask_hex") != hex(actual_mask):
                fail(errors, f"E={E} code={code}: hexadecimal bitmask mismatch")
            try:
                ap, bp, endpoint_dev = reconstruct(groups, selected, sorted(expected_ids))
                product = conv(ap, bp)
                if len(product) != len(target):
                    fail(errors, f"E={E} code={code}: product degree mismatch")
                    continue
                residual = max(abs(x - y) for x, y in zip(product, target))
                minimum, maximum, score, fractional = metric(ap, bp, eta)
                max_residual = max(max_residual, residual)
                max_endpoint_deviation = max(max_endpoint_deviation, endpoint_dev)
                if residual > tol:
                    fail(errors, f"E={E} code={code}: MP residual {mp.nstr(residual, 8)} > {args.tolerance}")
                if endpoint_dev > tol:
                    fail(errors, f"E={E} code={code}: pre-pin endpoint deviation {mp.nstr(endpoint_dev, 8)} > {args.tolerance}")
                if score < 0:
                    fail(errors, f"E={E} code={code}: negative violation score")
                if minimum >= 0 and maximum <= 1 and fractional:
                    closed_rows += 1
                saved = item.get("cpu_float64_independent_root_reconstruction")
                if saved is None:
                    missing_cpu_rows += 1
                else:
                    cpu_copy_rows += 1
                    for msg in cpu_copy_compare(saved, ap, bp, residual, minimum, maximum, fractional, tol):
                        fail(errors, f"E={E} code={code}: {msg}")
            except Exception as exc:  # Make any row-level reconstruction failure visible in receipt.
                fail(errors, f"E={E} code={code}: reconstruction exception {type(exc).__name__}: {exc}")
    complete = checkpoint.get("status") == "complete" and all(x["status"] == "complete" and x["coverage"][1] == x["total_codes"] for x in per_E)
    if not complete and not args.allow_partial:
        fail(errors, "checkpoint is incomplete; rerun with --allow-partial to audit only its saved prefix")
    status = "PASS_COMPLETE" if complete and not errors else ("PASS_PARTIAL_PREFIX" if not errors else "FAIL")
    result = {
        "audit": "independent-R8-four-block-v2-MP80",
        "result": status,
        "claim": ("saved-prefix consistency only; no exclusion conclusion" if not complete
                  else "ledger consistency only; no general mathematical conclusion"),
        "checkpoint_status": checkpoint.get("status"),
        "a": a, "b": b, "n_blocks": blocks,
        "franel_layout_total": expected_layout_total(a, b, blocks),
        "franel_formula": "F_4(a)F_4(b)/2 + 2^(4-1) - 2",
        "hashes": hashes,
        "mp_dps": args.mp_dps, "tolerance": str(args.tolerance),
        "rows_reconstructed": rows, "rows_with_stored_cpu_copy": cpu_copy_rows,
        "rows_missing_stored_cpu_copy": missing_cpu_rows,
        "closed_unit_interval_fractional_saved_rows": closed_rows,
        "max_mp_residual": mp.nstr(max_residual, 30),
        "max_endpoint_deviation_without_pinning": mp.nstr(max_endpoint_deviation, 30),
        "max_lift_root_to_float64_quadratic_delta": mp.nstr(max_lift_root_float64_delta, 30),
        "per_E": per_E,
        "python": sys.version, "platform": platform.platform(), "mpmath": mp.__version__,
        "errors": errors,
    }
    return result, (0 if not errors else 1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--options", required=True)
    parser.add_argument("--roots-input", required=True)
    parser.add_argument("--engine-source", required=True)
    parser.add_argument("--record", required=True)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--tolerance", type=float, default=1e-7)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    result, rc = audit(args)
    out = Path(args.record)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print(f"{result['result']}: {out}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
