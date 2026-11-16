#!/usr/bin/env python3
"""General streaming balanced-block MPS double-single numerical screen.

The prepared options contain a root-balanced split in each block.  This tool
uses Metal double-single convolution only to rank finite choices.  Root-factor
float64 reconstruction is deliberately independent and is retained for every
small smoke choice and the final top 64 choices per E.  Neither check is an
algebraic certification of the root data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import torch


METAL = r'''#pragma clang fp contract(off)
#include <metal_stdlib>
using namespace metal;
inline float2 qsum(float a,float b){float s=a+b;return float2(s,b-(s-a));}
inline float2 tsum(float a,float b){float s=a+b,bb=s-a;return float2(s,(a-(s-bb))+(b-bb));}
inline float2 tprod(float a,float b){float p=a*b;return float2(p,fma(a,b,-p));}
inline float2 dsadd(float2 a,float2 b){float2 s=tsum(a.x,b.x);return qsum(s.x,s.y+a.y+b.y);}
inline float2 dsmul(float2 a,float2 b){float2 p=tprod(a.x,b.x);return qsum(p.x,p.y+a.x*b.y+a.y*b.x+a.y*b.y);}
kernel void dsconv(device float* oh [[buffer(0)]], device float* ol [[buffer(1)]],
 device const float* ah [[buffer(2)]], device const float* al [[buffer(3)]],
 device const float* bh [[buffer(4)]], device const float* bl [[buffer(5)]],
 constant uint& aw [[buffer(6)]], constant uint& bw [[buffer(7)]], uint gid [[thread_position_in_grid]]){
 uint ow=aw+bw-1,row=gid/ow,k=gid%ow,lo=(k>=bw-1)?k-(bw-1):0,hi=min(k,aw-1);float2 z=float2(0.0f);
 for(uint i=lo;i<=hi;i++)z=dsadd(z,dsmul(float2(ah[row*aw+i],al[row*aw+i]),float2(bh[row*bw+k-i],bl[row*bw+k-i])));
 oh[gid]=z.x;ol[gid]=z.y;
}
'''


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def dsconv(lib, a, b):
    """One batched DS convolution; inputs are (hi, lo, coefficient width)."""
    ah, al, aw = a
    bh, bl, bw = b
    rows = ah.numel() // aw
    ow = aw + bw - 1
    oh = torch.empty(rows * ow, dtype=torch.float32, device="mps")
    ol = torch.empty_like(oh)
    lib.dsconv(oh, ol, ah, al, bh, bl, aw, bw)
    return oh, ol, ow


def selected_groups(blocks: list[dict], E: int, L: int, digits: list[int]) -> tuple[int, list[dict]]:
    """Recover an auditable root-group mask without trusting block products."""
    mask = 0
    rows = []
    for block, digit in zip(blocks, digits):
        option = block["options"][str(E)][digit]
        if E < L:
            picks = [block["lift_groups"][q] for q in option["selected_lift_positions"]]
            mode = "selected_lifts"
        elif option["mode"] == "original":
            picks, mode = [block["original_group"]], "original"
        else:
            picks, mode = list(block["lift_groups"]), "all_lifts"
        mask |= sum(1 << int(q) for q in picks)
        rows.append({"block": int(block["block"]), "option_index": int(digit), "mode": mode,
                     "A_group_indices": [int(q) for q in picks]})
    return mask, rows


def roots_reconstruct(groups: list[np.ndarray], mask: int) -> tuple[np.ndarray, np.ndarray]:
    """Independent float64 product of individual root groups."""
    a = np.array([1.0])
    b = np.array([1.0])
    for j, factor in enumerate(groups):
        if (mask >> j) & 1:
            a = np.convolve(a, factor)
        else:
            b = np.convolve(b, factor)
    a /= a[0]
    b /= b[0]
    # The option generator's root-balance identity fixes these endpoints.
    a[0] = a[-1] = 1.0
    b[0] = b[-1] = 1.0
    return a, b


def gpu_poly(h: np.ndarray, l: np.ndarray, row: int, width: int, degree: int) -> np.ndarray:
    values = h[row * width:(row + 1) * width].astype(np.float64)
    values += l[row * width:(row + 1) * width].astype(np.float64)
    values = values[:degree + 1] / values[0]
    values[0] = values[-1] = 1.0
    return values


def choose_options(data: dict, E: int, L: int) -> list[list[dict]]:
    options = [list(block["options"][str(E)]) for block in data["blocks"]]
    if E == L:
        # Complement representatives: the first parent is fixed to original.
        options[0] = [options[0][0]]
    return options


def option_layouts(data: dict, L: int) -> dict[int, tuple[list[int], int]]:
    """Return the canonical mixed-radix layout for every available E."""
    if not data["blocks"]:
        raise ValueError("options file has no blocks")
    keys = set(data["blocks"][0]["options"])
    if not keys:
        raise ValueError("options file has no E choices")
    for block in data["blocks"][1:]:
        if set(block["options"]) != keys:
            raise ValueError("blocks do not provide the same E choices")
    layouts = {}
    for text_E in keys:
        E = int(text_E)
        if str(E) != text_E:
            raise ValueError(f"non-canonical E key {text_E!r}")
        radices = [len(choices) for choices in choose_options(data, E, L)]
        if any(radix <= 0 for radix in radices):
            raise ValueError(f"E={E} has an empty option list")
        layouts[E] = (radices, math.prod(radices))
    return layouts


def validate_resume_entries(entries: list[dict], layouts: dict[int, tuple[list[int], int]]) -> dict[int, dict]:
    """Validate every persisted ledger row before any checkpoint write."""
    prior: dict[int, dict] = {}
    for entry in entries:
        E = int(entry.get("E"))
        if E in prior or E not in layouts:
            raise ValueError(f"checkpoint has duplicate or unavailable E={E}")
        radices, total = layouts[E]
        if entry.get("radices_block0_first") != radices or entry.get("total_codes") != total:
            raise ValueError(f"checkpoint E={E} mixed-radix layout differs from current input")
        coverage = entry.get("coverage_completed_half_open")
        if (not isinstance(coverage, list) or len(coverage) != 2 or
                type(coverage[0]) is not int or type(coverage[1]) is not int or
                coverage[0] != 0 or not 0 <= coverage[1] <= total):
            raise ValueError(f"checkpoint E={E} has invalid contiguous coverage")
        counts = entry.get("counts")
        if not isinstance(counts, dict) or type(counts.get("codes_seen")) is not int:
            raise ValueError(f"checkpoint E={E} has invalid counts")
        if counts["codes_seen"] != coverage[1]:
            raise ValueError(f"checkpoint E={E} codes_seen disagrees with coverage")
        for key in ("fractional_masks", "closed_unit_interval_fractional_masks"):
            if type(counts.get(key)) is not int or not 0 <= counts[key] <= counts["codes_seen"]:
                raise ValueError(f"checkpoint E={E} has invalid {key}")
        if counts["closed_unit_interval_fractional_masks"] > counts["fractional_masks"]:
            raise ValueError(f"checkpoint E={E} has incompatible mask counts")
        prior[E] = entry
    return prior


def resume_start(entry: dict, requested_end: int) -> int:
    """Allow coverage extension or an equal bound, never ledger shrinkage."""
    covered = int(entry.get("coverage_completed_half_open", [0, 0])[1]) if entry else 0
    if covered > requested_end:
        raise ValueError(f"requested end {requested_end} shrinks existing coverage {covered}")
    return covered


def prepack_options(options: list[list[dict]]) -> dict:
    """Pad every option once, then transfer the option tables to MPS once.

    Batch processing uses only GPU index gathers from these tables: it never
    runs a Python loop over a row and never calls np.pad in its hot path.
    """
    packed: dict[str, list[dict]] = {"Ablock": [], "Bblock": []}
    for side in packed:
        for block_options in options:
            width = max(len(choice[side]["hi"]) for choice in block_options)
            count = len(block_options)
            hi = np.zeros((count, width), dtype=np.float32)
            lo = np.zeros_like(hi)
            degrees = np.empty(count, dtype=np.int64)
            for i, choice in enumerate(block_options):
                h = np.asarray(choice[side]["hi"], dtype=np.float32)
                low = np.asarray(choice[side]["lo"], dtype=np.float32)
                if h.shape != low.shape:
                    raise ValueError(f"mismatched hi/lo shapes in {side}")
                hi[i, :h.size] = h
                lo[i, :low.size] = low
                degrees[i] = h.size - 1
            packed[side].append({
                "hi": torch.from_numpy(hi).to("mps"),
                "lo": torch.from_numpy(lo).to("mps"),
                "degree": torch.from_numpy(degrees).to("mps"),
                "width": width,
            })
    return packed


def digits_for_codes(codes: np.ndarray, radices: list[int]) -> np.ndarray:
    q = codes.copy()
    digits = np.empty((codes.size, len(radices)), dtype=np.int64)
    for j, radix in enumerate(radices):
        digits[:, j] = q % radix
        q //= radix
    return digits


def gathered_product(lib, tables: list[dict], device_digits: list[torch.Tensor]):
    gathered = []
    degrees = []
    for table, index in zip(tables, device_digits):
        # index_select is vectorized MPS gather over a prepacked option table.
        gathered.append((table["hi"].index_select(0, index).reshape(-1),
                         table["lo"].index_select(0, index).reshape(-1), table["width"]))
        degrees.append(table["degree"].index_select(0, index))
    product = gathered[0]
    for other in gathered[1:]:
        product = dsconv(lib, product, other)
    return product, torch.stack(degrees).sum(dim=0)


def run_batch(lib, packed: dict, codes: np.ndarray, radices: list[int], eta: float,
              materialize_coefficients: bool = False) -> dict:
    """Run 2*(number_of_blocks-1) Metal convolutions and copy only metrics.

    Coefficients cross back to CPU exclusively for smoke reconstruction.
    """
    digits = digits_for_codes(codes, radices)
    device_digits = [torch.from_numpy(digits[:, j]).to("mps") for j in range(digits.shape[1])]
    (ah, al, aw), da = gathered_product(lib, packed["Ablock"], device_digits)
    (bh, bl, bw), db = gathered_product(lib, packed["Bblock"], device_digits)
    rows = codes.size
    av = (ah + al).reshape(rows, aw)
    bv = (bh + bl).reshape(rows, bw)
    av = av / av[:, :1]
    bv = bv / bv[:, :1]
    ia = torch.arange(aw, device="mps")[None, :]
    ib = torch.arange(bw, device="mps")[None, :]
    interior_a = (ia > 0) & (ia < da[:, None])
    interior_b = (ib > 0) & (ib < db[:, None])
    # Explicitly restore normalized endpoints; padded tails are excluded.
    av = torch.where((ia == 0) | (ia == da[:, None]), torch.ones_like(av), av)
    bv = torch.where((ib == 0) | (ib == db[:, None]), torch.ones_like(bv), bv)
    active = torch.cat((interior_a, interior_b), dim=1)
    values = torch.cat((av, bv), dim=1)
    score = (torch.relu(-values) + torch.relu(values - 1.0)).masked_fill(~active, 0.0).sum(dim=1)
    guarded = torch.where(active, values, torch.ones_like(values))
    mn, mx = guarded.min(dim=1).values, guarded.max(dim=1).values
    maxabs = torch.where(active, values.abs(), torch.ones_like(values)).max(dim=1).values
    fractional = ((values > eta) & (values < 1.0 - eta) & active).any(dim=1)
    torch.mps.synchronize()
    result = {
        "digits": digits,
        "degree_A": da.cpu().numpy(), "degree_B": db.cpu().numpy(),
        "width_A": aw, "width_B": bw,
        "score": score.cpu().numpy(), "min": mn.cpu().numpy(), "max": mx.cpu().numpy(),
        "maxabs": maxabs.cpu().numpy(), "fractional": fractional.cpu().numpy().astype(bool),
    }
    for name in ("score", "min", "max", "maxabs"):
        if not np.all(np.isfinite(result[name])):
            raise FloatingPointError(f"non-finite MPS metric {name}; refusing to update checkpoint counts")
    if materialize_coefficients:
        result.update({"Ahi": ah.cpu().numpy(), "Alo": al.cpu().numpy(),
                       "Bhi": bh.cpu().numpy(), "Blo": bl.cpu().numpy()})
    del ah, al, bh, bl, da, db, av, bv, ia, ib, interior_a, interior_b
    del active, values, score, guarded, mn, mx, maxabs, fractional, device_digits
    return result


def candidate_record(data: dict, E: int, L: int, code: int, digits: list[int], gpu: dict, row: int) -> dict:
    mask, mapping = selected_groups(data["blocks"], E, L, digits)
    return {
        "code": int(code), "mixed_radix_digits_block0_first": [int(q) for q in digits],
        "A_root_group_bitmask_decimal": str(mask), "A_root_group_bitmask_hex": hex(mask),
        "mask_mapping": mapping, "degree_A": int(gpu["degree_A"][row]), "degree_B": int(gpu["degree_B"][row]),
        "gpu_ds_violation_l1": float(gpu["score"][row]),
        "gpu_ds_min_interior": float(gpu["min"][row]), "gpu_ds_max_interior": float(gpu["max"][row]),
        "gpu_ds_max_abs_interior": float(gpu["maxabs"][row]),
        "gpu_ds_has_fractionality": bool(gpu["fractional"][row]),
    }


def recheck_top(entry: dict, groups: list[np.ndarray], C: np.ndarray) -> None:
    target = C / C[0]
    for row in entry["top64"]:
        a, b = roots_reconstruct(groups, int(row["A_root_group_bitmask_decimal"]))
        row["cpu_float64_independent_root_reconstruction"] = {
            "A_coefficients_ascending": a.tolist(), "B_coefficients_ascending": b.tolist(),
            "convolution_max_abs_residual": float(np.max(np.abs(np.convolve(a, b) - target))),
            "minimum_actual_coefficient": float(min(a.min(), b.min())),
            "maximum_actual_coefficient": float(max(a.max(), b.max())),
            "has_fractionality_eta_1e-5": bool(np.any((a[1:-1] > 1e-5) & (a[1:-1] < 1 - 1e-5)) or
                                                 np.any((b[1:-1] > 1e-5) & (b[1:-1] < 1 - 1e-5))),
        }


def exact_fair_check(raw: dict, gpu: dict, row: int, L: int) -> dict:
    h = np.asarray(raw["H_coefficients_ascending"], dtype=np.float64)
    h_lift = np.zeros((len(h) - 1) * L + 1)
    h_lift[::L] = h
    ga = gpu_poly(gpu["Ahi"], gpu["Alo"], row, gpu["width_A"], int(gpu["degree_A"][row]))
    gb = gpu_poly(gpu["Bhi"], gpu["Blo"], row, gpu["width_B"], int(gpu["degree_B"][row]))
    return {"code": 0, "A_vs_exact_integer_H_max_abs": float(np.max(np.abs(ga - h))),
            "B_vs_exact_integer_H_of_xL_max_abs": float(np.max(np.abs(gb - h_lift)))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--options", type=Path, required=True)
    parser.add_argument("--roots-input", type=Path, required=True)
    parser.add_argument("--candidate-index", type=int, default=0)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=32768)
    parser.add_argument("--max-codes", type=int)
    parser.add_argument("--smoke", action="store_true", help="run E=1,4,L at at most 1024 codes each")
    parser.add_argument("--benchmark-one-batch", action="store_true", help="time exactly one 32768-code E=4 batch")
    parser.add_argument("--eta", type=float, default=1e-5)
    args = parser.parse_args()
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data is not mounted; refusing to write elsewhere")
    if args.batch_size <= 0 or args.batch_size > 131072:
        raise SystemExit("batch size must be 1..131072")
    if not torch.backends.mps.is_available():
        raise SystemExit("MPS unavailable")
    data = json.loads(args.options.read_text())
    raw = json.loads(args.roots_input.read_text())["candidates"][args.candidate_index]
    L = int(data["L"])
    block_count = int(data.get("block_count", len(data["blocks"])))
    if (data.get("candidate") != raw.get("name") or L != int(raw["L"]) or
            len(data["blocks"]) != block_count or int(raw.get("parent_block_count", block_count)) != block_count):
        raise SystemExit("options/roots candidate, L, or block-count mismatch")
    if block_count < 2:
        raise SystemExit("need at least two blocks for streaming DS convolution")
    try:
        layouts = option_layouts(data, L)
    except ValueError as exc:
        raise SystemExit(f"invalid options layout: {exc}") from exc
    groups = [np.asarray(g["coefficients_ascending"], dtype=np.float64) for g in raw["groups"]]
    C = np.asarray(raw["C_coefficients_ascending"], dtype=np.float64)
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "format": "streaming-balanced-double-single-mps-v2-general-blocks",
        "status": "running", "source_sha256": sha(Path(__file__)), "options_sha256": sha(args.options),
        "roots_input_sha256": sha(args.roots_input), "candidate": data["candidate"],
        "candidate_index": args.candidate_index,
        "parameters": {"batch_size": args.batch_size, "eta": args.eta, "max_codes": args.max_codes,
                       "smoke": args.smoke, "benchmark_one_batch": args.benchmark_one_batch},
        "algorithm": {"device": "mps", "float_format": "double-single float32 hi+lo", "L": L,
                      "local_block_count": block_count, "batch_DS_convolutions": 2 * (block_count - 1),
                      "mixed_radix": "block 0 least significant", "endpoint_normalization": "endpoints explicitly set to 1; padded tails excluded",
                      "hot_path": "prepacked MPS option tables plus vectorized index_select gathers; metrics/top64 copied to CPU"},
        "input_identity": {"H_coefficients_ascending": raw.get("H_coefficients_ascending"),
                           "balanced_enumeration": raw.get("balanced_enumeration"),
                           "all_balanced_nonconstant_complement_representatives": raw.get("balanced_enumeration", {}).get("nonconstant_complement_representatives")},
        "versions": {"python": sys.version, "torch": torch.__version__, "platform": platform.platform()},
        "E_progress": [], "claim_boundary": "MPS double-single ranks finite numerical candidates; it is not a proof of acceptance or exclusion.",
    }
    selected_es = [E for E in (1, 4, L) if E in layouts] if args.smoke else sorted(layouts)

    def requested_end(E: int) -> int:
        total = layouts[E][1]
        limit = 1024 if args.smoke else (args.max_codes if args.max_codes is not None else total)
        return min(total, limit)

    prior: dict[int, dict] = {}
    if args.checkpoint.exists():
        old = json.loads(args.checkpoint.read_text())
        keys = ("format", "source_sha256", "options_sha256", "roots_input_sha256", "candidate", "candidate_index")
        old_parameters = old.get("parameters", {})
        # eta changes the fractional and unit-interval counts, hence it is a
        # scoring/calibration identity rather than a tunable runtime setting.
        # batch_size and max_codes intentionally remain resumable: they only
        # alter chunking and requested contiguous coverage.
        if (not all(old.get(key) == report.get(key) for key in keys) or
                old_parameters.get("smoke") != args.smoke or
                old_parameters.get("eta") != args.eta):
            raise SystemExit("checkpoint exists but source/input/smoke/eta identity differs; refusing to overwrite it")
        try:
            prior = validate_resume_entries(old.get("E_progress", []), layouts)
        except (TypeError, ValueError) as exc:
            raise SystemExit(f"checkpoint ledger validation failed: {exc}") from exc
        # Seed every old E before a benchmark or an earlier E can atomically
        # write.  Replacing only the current E below cannot erase later work.
        report["E_progress"] = list(old.get("E_progress", []))
        if "benchmark_one_32768_code_batch" in old:
            report["benchmark_one_32768_code_batch"] = old["benchmark_one_32768_code_batch"]
        report["resumed_from_checkpoint"] = str(args.checkpoint)
        # Check all requested rows before a benchmark checkpoint write.  An
        # accidental lower --max-codes must fail rather than partly run and
        # only later discover a ledger-shrinking request.
        try:
            for E in selected_es:
                if E in prior:
                    resume_start(prior[E], requested_end(E))
        except ValueError as exc:
            raise SystemExit(f"checkpoint resume bound validation failed: {exc}") from exc
    started = perf_counter()
    try:
        lib = torch.mps.compile_shader(METAL)
        prepared: dict[int, tuple[list[list[dict]], list[int], dict]] = {}
        def get_prepared(E: int):
            if E not in prepared:
                opts = choose_options(data, E, L)
                radices = [len(x) for x in opts]
                prepared[E] = (opts, radices, prepack_options(opts))
            return prepared[E]
        if args.benchmark_one_batch:
            if "4" not in data["blocks"][0]["options"]:
                raise SystemExit("E=4 is unavailable for the requested benchmark")
            _, radices, packed = get_prepared(4)
            if math.prod(radices) < 32768:
                raise SystemExit("E=4 has fewer than 32768 codes; requested one-batch benchmark is undefined")
            t0 = perf_counter()
            bench = run_batch(lib, packed, np.arange(32768, dtype=np.int64), radices, args.eta)
            elapsed = perf_counter() - t0
            report["benchmark_one_32768_code_batch"] = {"E": 4, "batch_size": 32768, "elapsed_seconds": elapsed,
                                                         "codes_per_second": 32768 / elapsed,
                                                         "coefficient_materialized_to_cpu": False,
                                                         "metric_score_min": float(bench["score"].min())}
            del bench
            torch.mps.empty_cache()
            atomic(args.checkpoint, report)
        es = selected_es
        for E in dict.fromkeys(es):
            _, radices, packed = get_prepared(E)
            total = math.prod(radices)
            end = requested_end(E)
            old = prior.get(E, {})
            try:
                start = resume_start(old, end)
            except ValueError as exc:
                raise SystemExit(f"checkpoint E={E}: {exc}") from exc
            best = list(old.get("top64", []))
            counts = dict(old.get("counts", {"codes_seen": 0, "fractional_masks": 0, "closed_unit_interval_fractional_masks": 0}))
            previous_smoke = old.get("smoke_gpu_vs_independent_cpu_float64_root_reconstruction", {})
            smoke_errors: list[float] = list(previous_smoke.get("per_code_errors", [])) if args.smoke else []
            fair = old.get("fair_EL_exact_integer_coefficient_check") if args.smoke else None
            for lo in range(start, end, args.batch_size):
                hi = min(end, lo + args.batch_size)
                codes = np.arange(lo, hi, dtype=np.int64)
                gpu = run_batch(lib, packed, codes, radices, args.eta, materialize_coefficients=args.smoke)
                counts["codes_seen"] += len(codes)
                counts["fractional_masks"] += int(gpu["fractional"].sum())
                counts["closed_unit_interval_fractional_masks"] += int((gpu["fractional"] & (gpu["min"] >= 0) & (gpu["max"] <= 1)).sum())
                pick = np.lexsort((codes, gpu["score"]))[:min(64, len(codes))]
                best.extend(candidate_record(data, E, L, int(codes[i]), gpu["digits"][i].tolist(), gpu, int(i)) for i in pick)
                best = sorted(best, key=lambda row: (row["gpu_ds_violation_l1"], row["code"]))[:64]
                if args.smoke:
                    for i, code in enumerate(codes):
                        mask, _ = selected_groups(data["blocks"], E, L, gpu["digits"][i].tolist())
                        ca, cb = roots_reconstruct(groups, mask)
                        ga = gpu_poly(gpu["Ahi"], gpu["Alo"], i, gpu["width_A"], int(gpu["degree_A"][i]))
                        gb = gpu_poly(gpu["Bhi"], gpu["Blo"], i, gpu["width_B"], int(gpu["degree_B"][i]))
                        smoke_errors.append(max(float(np.max(np.abs(ga - ca))), float(np.max(np.abs(gb - cb)))))
                        if E == L and int(code) == 0:
                            fair = exact_fair_check(raw, gpu, i, L)
                entry = {"E": E, "status": "running", "radices_block0_first": radices, "total_codes": total,
                         "coverage_completed_half_open": [0, hi], "counts": counts, "top64": best}
                report["E_progress"] = [row for row in report["E_progress"] if row["E"] != E] + [entry]
                report["elapsed_seconds"] = perf_counter() - started
                atomic(args.checkpoint, report)
                del gpu
                torch.mps.empty_cache()
            entry = {"E": E, "status": "complete" if end == total else "partial_unknown", "radices_block0_first": radices,
                     "total_codes": total, "coverage_completed_half_open": [0, end], "counts": counts, "top64": best}
            if args.smoke:
                entry["smoke_gpu_vs_independent_cpu_float64_root_reconstruction"] = {"codes_checked": end,
                    "max_abs_error": float(max(smoke_errors, default=0.0)), "per_code_errors": smoke_errors}
                if fair is not None:
                    entry["fair_EL_exact_integer_coefficient_check"] = fair
            recheck_top(entry, groups, C)
            report["E_progress"] = [row for row in report["E_progress"] if row["E"] != E] + [entry]
            report["elapsed_seconds"] = perf_counter() - started
            atomic(args.checkpoint, report)
        report["status"] = "smoke_complete" if args.smoke else ("complete" if all(row["status"] == "complete" for row in report["E_progress"]) else "partial_unknown")
        report["elapsed_seconds"] = perf_counter() - started
        atomic(args.checkpoint, report)
        print(json.dumps({"status": report["status"], "checkpoint": str(args.checkpoint), "E": es}, indent=2))
    finally:
        if torch.backends.mps.is_available():
            torch.mps.synchronize()
            torch.mps.empty_cache()


if __name__ == "__main__":
    main()
