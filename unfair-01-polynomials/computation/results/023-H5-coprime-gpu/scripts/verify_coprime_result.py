#!/usr/bin/env python3
"""Independent CPU verifier for saved coprime DS search results.

This program never invokes MPS or imports the search engine.  It rebuilds
mixed-radix choices from the option cache's explicit selected_group_indices,
then recomputes every retained top64 root product with mpmath at configurable
precision.  It is numerical evidence, not a root-isolation certificate.
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


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def franel(n: int) -> int:
    return sum(math.comb(n, k) ** 3 for k in range(n + 1))


def dilate(h: list[int], stride: int) -> list[int]:
    out = [0] * ((len(h) - 1) * stride + 1)
    out[::stride] = h
    return out


def integer_convolve(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def choose_options(data: dict, E: int, sub_a: int, sub_b: int) -> list[list[dict]]:
    choices = [list(block["options"][str(E)]) for block in data["blocks"]]
    if E == sub_a * sub_b:
        fixed = [entry for entry in choices[0] if entry.get("mode") == "all_a"]
        if len(fixed) != 1:
            raise ValueError("terminal E has no unique block-0 all_a option")
        choices[0] = fixed
    return choices


def layouts(data: dict, sub_a: int, sub_b: int) -> dict[int, tuple[list[int], int]]:
    keyset = set(data["blocks"][0]["options"])
    if not keyset:
        raise ValueError("empty option cache")
    for block in data["blocks"][1:]:
        if set(block["options"]) != keyset:
            raise ValueError("E keys differ between blocks")
    answer = {}
    for text in keyset:
        E = int(text)
        if str(E) != text or not 0 < E <= sub_a * sub_b:
            raise ValueError(f"invalid E key {text}")
        radices = [len(x) for x in choose_options(data, E, sub_a, sub_b)]
        if any(r <= 0 for r in radices):
            raise ValueError(f"empty options at E={E}")
        answer[E] = (radices, math.prod(radices))
    expected = franel(sub_a) * franel(sub_b) // 2 + 2
    if sum(total for _, total in answer.values()) != expected:
        raise ValueError("layout total does not match Franel formula")
    return answer


def decode(code: int, radices: list[int]) -> list[int]:
    if not 0 <= code < math.prod(radices):
        raise ValueError("code outside mixed-radix range")
    result = []
    for radix in radices:
        result.append(code % radix)
        code //= radix
    return result


def expected_mapping(data: dict, E: int, digits: list[int]) -> tuple[set[int], list[dict]]:
    sub_a, sub_b = int(data["a"]), int(data["b"])
    options = choose_options(data, E, sub_a, sub_b)
    selected: set[int] = set()
    mapping = []
    for block, block_options, digit in zip(data["blocks"], options, digits):
        option = block_options[digit]
        picked = [int(x) for x in option["selected_group_indices"]]
        if len(picked) != len(set(picked)):
            raise ValueError("cached option repeats a selected group")
        selected.update(picked)
        mapping.append({"block": int(block["block"]), "option_index": int(digit),
                        "mode": str(option["mode"]), "A_group_indices": picked})
    return selected, mapping


def poly_mul(a: list[mp.mpf], b: list[mp.mpf]) -> list[mp.mpf]:
    out = [mp.mpf("0") for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def reconstruct(groups: list[dict], selected: set[int]) -> tuple[list[mp.mpf], list[mp.mpf], tuple[mp.mpf, mp.mpf]]:
    a, b = [mp.mpf("1")], [mp.mpf("1")]
    for index, group in enumerate(groups):
        factor = [mp.mpf(str(value)) for value in group["coefficients_ascending"]]
        if index in selected:
            a = poly_mul(a, factor)
        else:
            b = poly_mul(b, factor)
    a = [x / a[0] for x in a]
    b = [x / b[0] for x in b]
    endpoint_deviation = (abs(a[-1] - 1), abs(b[-1] - 1))
    # The screen defines the exact balanced endpoints as one after normalization.
    a[0] = a[-1] = mp.mpf("1")
    b[0] = b[-1] = mp.mpf("1")
    return a, b, endpoint_deviation


def mp_metrics(a: list[mp.mpf], b: list[mp.mpf], eta: mp.mpf) -> dict[str, object]:
    interior = a[1:-1] + b[1:-1]
    score = sum(max(-x, mp.mpf("0")) + max(x - 1, mp.mpf("0")) for x in interior)
    return {"score": score, "minimum": min(interior), "maximum": max(interior),
            "fractional": any(eta < x < 1 - eta for x in interior)}


def max_difference(mp_values: list[mp.mpf], saved: list[float]) -> mp.mpf:
    if len(mp_values) != len(saved):
        raise AssertionError("saved polynomial length differs from independent reconstruction")
    return max(abs(x - mp.mpf(str(y))) for x, y in zip(mp_values, saved))


def validate_cache_options(data: dict, raw: dict) -> None:
    sub_a, sub_b = int(data["a"]), int(data["b"])
    all_groups = set(range(len(raw["groups"])))
    for block in data["blocks"]:
        a_groups = [int(x) for x in block["a_lift_groups"]]
        b_groups = [int(x) for x in block["b_lift_groups"]]
        if len(a_groups) != sub_a or len(b_groups) != sub_b or set(a_groups) & set(b_groups):
            raise ValueError("bad per-block a/b lift grouping")
        if set(a_groups + b_groups) != set(map(int, block["all_group_indices"])):
            raise ValueError("per-block group union mismatch")
        for text_E, entries in block["options"].items():
            E = int(text_E)
            for option in entries:
                picked = [int(x) for x in option["selected_group_indices"]]
                if len(picked) != len(set(picked)) or not set(picked) <= set(a_groups + b_groups):
                    raise ValueError("option selected groups are invalid")
                if int(option["u"]) * sub_b + int(option["v"]) * sub_a != E:
                    raise ValueError("option violates E=b*u+a*v")
                if option["mode"] == "all_a" and picked != a_groups:
                    raise ValueError("all_a mapping is not the exact a-lift group list")
                if option["mode"] == "all_b" and picked != b_groups:
                    raise ValueError("all_b mapping is not the exact b-lift group list")
                for side in ("Ablock", "Bblock"):
                    part = option[side]
                    if not (len(part["hi"]) == len(part["lo"]) == len(part["float64"])):
                        raise ValueError("option hi/lo/float64 widths differ")
    if set(range(len(raw["groups"]))) != all_groups:
        raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--options", type=Path, required=True)
    parser.add_argument("--roots-input", type=Path, required=True)
    parser.add_argument("--engine-source", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--tolerance", type=str, default="1e-7")
    parser.add_argument("--allow-partial", action="store_true", help="audit stored contiguous prefixes, but do not claim full coverage")
    args = parser.parse_args()
    if args.mp_dps < 50:
        raise SystemExit("mp-dps must be at least 50")
    mp.mp.dps = args.mp_dps
    tolerance = mp.mpf(args.tolerance)
    checkpoint = json.loads(args.checkpoint.read_text())
    data = json.loads(args.options.read_text())
    raw = json.loads(args.roots_input.read_text())["candidates"][int(checkpoint["candidate_index"])]
    if checkpoint.get("format") != "coprime-substitution-double-single-mps-v1":
        raise SystemExit("checkpoint is not from the frozen coprime engine format")
    if (checkpoint.get("source_sha256") != sha(args.engine_source) or checkpoint.get("options_sha256") != sha(args.options) or
            checkpoint.get("roots_input_sha256") != sha(args.roots_input) or checkpoint.get("candidate") != raw["name"]):
        raise SystemExit("source/options/roots hash or candidate binding mismatch")
    if raw.get("kind") != "coprime_substitution" or data.get("kind") != "coprime_substitution":
        raise SystemExit("input/options kind mismatch")
    sub_a, sub_b = int(raw["a"]), int(raw["b"])
    if (int(data["a"]), int(data["b"])) != (sub_a, sub_b) or math.gcd(sub_a, sub_b) != 1:
        raise SystemExit("coprime a/b mismatch")
    if int(raw["balance_terminal_E"]) != sub_a * sub_b or int(data["balance_terminal_E"]) != sub_a * sub_b:
        raise SystemExit("terminal E mismatch")
    if len(data["blocks"]) != len(raw["blocks"]) != 3:
        raise SystemExit("expected exactly three compatible blocks")
    h = [int(x) for x in raw["H_coefficients_ascending"]]
    if raw["C_coefficients_ascending"] != integer_convolve(dilate(h, sub_a), dilate(h, sub_b)):
        raise SystemExit("raw C is not exact H(x^a)H(x^b)")
    validate_cache_options(data, raw)
    all_layouts = layouts(data, sub_a, sub_b)
    progress = checkpoint.get("E_progress", [])
    progress_by_E = {int(entry["E"]): entry for entry in progress}
    if len(progress_by_E) != len(progress):
        raise SystemExit("duplicate E progress entries")
    full = checkpoint.get("status") == "complete" and set(progress_by_E) == set(all_layouts)
    if not full and not args.allow_partial:
        raise SystemExit("checkpoint is partial/non-complete; rerun only with --allow-partial")
    eta = mp.mpf(str(checkpoint["parameters"]["eta"]))
    audited_rows = 0
    max_saved_poly_difference = mp.mpf("0")
    max_conv_residual = mp.mpf("0")
    max_gpu_score_difference = mp.mpf("0")
    max_endpoint_deviation = mp.mpf("0")
    fair_checked = False
    e_summary = []
    for E, entry in sorted(progress_by_E.items()):
        if E not in all_layouts:
            raise SystemExit(f"checkpoint stores unavailable E={E}")
        radices, total = all_layouts[E]
        coverage = entry.get("coverage_completed_half_open")
        counts = entry.get("counts", {})
        if (coverage is None or coverage[0] != 0 or not isinstance(coverage[1], int) or
                not 0 <= coverage[1] <= total or entry.get("radices_block0_first") != radices or
                entry.get("total_codes") != total or counts.get("codes_seen") != coverage[1]):
            raise SystemExit(f"invalid coverage ledger for E={E}")
        if full and (coverage[1] != total or entry.get("status") != "complete"):
            raise SystemExit(f"full checkpoint has incomplete E={E}")
        rows = entry.get("top64", [])
        if len(rows) > 64:
            raise SystemExit(f"more than top64 saved for E={E}")
        scores = [float(row["gpu_ds_violation_l1"]) for row in rows]
        if scores != sorted(scores):
            raise SystemExit(f"top64 scores not sorted for E={E}")
        for row in rows:
            code = int(row["code"])
            if not 0 <= code < coverage[1]:
                raise SystemExit(f"saved code outside completed prefix at E={E}")
            digits = decode(code, radices)
            if [int(x) for x in row["mixed_radix_digits_block0_first"]] != digits:
                raise SystemExit(f"mixed-radix decode mismatch at E={E}, code={code}")
            selected, mapping = expected_mapping(data, E, digits)
            if row["mask_mapping"] != mapping:
                raise SystemExit(f"explicit selected_group_indices mapping mismatch at E={E}, code={code}")
            bitmask = sum(1 << index for index in selected)
            if row["A_root_group_bitmask_decimal"] != str(bitmask) or row["A_root_group_bitmask_hex"] != hex(bitmask):
                raise SystemExit(f"mask encoding mismatch at E={E}, code={code}")
            a_poly, b_poly, endpoint_dev = reconstruct(raw["groups"], selected)
            max_endpoint_deviation = max(max_endpoint_deviation, *endpoint_dev)
            saved = row.get("cpu_float64_independent_root_reconstruction")
            if not saved:
                raise SystemExit(f"top64 row lacks independent reconstruction at E={E}, code={code}")
            poly_diff = max(max_difference(a_poly, saved["A_coefficients_ascending"]),
                            max_difference(b_poly, saved["B_coefficients_ascending"]))
            if poly_diff > tolerance:
                raise SystemExit(f"MP{args.mp_dps} polynomial mismatch {poly_diff} exceeds tolerance at E={E}, code={code}")
            product = poly_mul(a_poly, b_poly)
            target = [mp.mpf(x) for x in raw["C_coefficients_ascending"]]
            residual = max(abs(x - y) for x, y in zip(product, target))
            if residual > tolerance:
                raise SystemExit(f"independent convolution residual {residual} exceeds tolerance at E={E}, code={code}")
            metrics = mp_metrics(a_poly, b_poly, eta)
            if bool(saved["has_fractionality_eta_1e-5"]) != bool(metrics["fractional"]):
                raise SystemExit(f"fractionality mismatch at E={E}, code={code}")
            max_saved_poly_difference = max(max_saved_poly_difference, poly_diff)
            max_conv_residual = max(max_conv_residual, residual)
            max_gpu_score_difference = max(max_gpu_score_difference, abs(metrics["score"] - mp.mpf(str(row["gpu_ds_violation_l1"]))))
            audited_rows += 1
            if E == sub_a * sub_b and code == 0:
                expected_a, expected_b = [mp.mpf(x) for x in dilate(h, sub_a)], [mp.mpf(x) for x in dilate(h, sub_b)]
                if max_difference(a_poly, expected_a) > tolerance or max_difference(b_poly, expected_b) > tolerance:
                    raise SystemExit("terminal code 0 does not reconstruct H(x^a), H(x^b)")
                if any(item["mode"] != "all_a" for item in mapping):
                    raise SystemExit("terminal code 0 is not all_a in every block")
                fair_checked = True
        e_summary.append({"E": E, "coverage": coverage, "total_codes": total, "top64_rows_checked": len(rows)})
    if full and not fair_checked:
        raise SystemExit("complete checkpoint lacks terminal fair code 0 in saved top64")
    record = {
        "status": "verified_complete" if full else "verified_partial_prefix_only",
        "claim_boundary": "MP reconstruction of saved rows is numerical evidence only; full status additionally validates only ledger coverage and not GPU arithmetic over unsaved codes.",
        "checkpoint": str(args.checkpoint), "checkpoint_sha256": sha(args.checkpoint),
        "engine_source": str(args.engine_source), "engine_source_sha256": sha(args.engine_source),
        "options_sha256": sha(args.options), "roots_input_sha256": sha(args.roots_input),
        "a": sub_a, "b": sub_b, "canonical_total": sum(total for _, total in all_layouts.values()),
        "e_summary": e_summary, "mp_dps": args.mp_dps, "tolerance": str(tolerance),
        "top64_rows_checked": audited_rows, "terminal_fair_code0_checked": fair_checked,
        "max_saved_polynomial_difference": str(max_saved_poly_difference),
        "max_independent_convolution_residual": str(max_conv_residual),
        "max_gpu_vs_mp_score_difference": str(max_gpu_score_difference),
        "max_raw_normalized_endpoint_deviation": str(max_endpoint_deviation),
        "python": sys.version, "platform": platform.platform(), "gpu_used": False,
    }
    args.record.parent.mkdir(parents=True, exist_ok=True)
    atomic(args.record, record)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
