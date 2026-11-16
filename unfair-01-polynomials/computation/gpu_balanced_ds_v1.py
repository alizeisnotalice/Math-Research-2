#!/usr/bin/env python3
"""Stream the L=12 balanced block choices with four MPS double-single convolutions.

This is a numerical screen, not an exclusion proof.  The input block options are
prepared at 100 decimal digits; the independent rechecks below deliberately go
back to the root-group input rather than multiplying those prepared options.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
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
inline float2 dsdiv(float2 a,float2 b){float q=a.x/b.x;float2 r=dsadd(a,-dsmul(float2(q,0.0f),b));return dsadd(float2(q,0.0f),float2((r.x+r.y)/b.x,0.0f));}
kernel void dsconv(device float* oh [[buffer(0)]], device float* ol [[buffer(1)]],
 device const float* ah [[buffer(2)]], device const float* al [[buffer(3)]],
 device const float* bh [[buffer(4)]], device const float* bl [[buffer(5)]],
 constant uint& aw [[buffer(6)]], constant uint& bw [[buffer(7)]], uint gid [[thread_position_in_grid]]){
 uint ow=aw+bw-1,row=gid/ow,k=gid%ow,lo=(k>=bw-1)?k-(bw-1):0,hi=min(k,aw-1);float2 z=float2(0.0f);
 for(uint i=lo;i<=hi;i++)z=dsadd(z,dsmul(float2(ah[row*aw+i],al[row*aw+i]),float2(bh[row*bw+k-i],bl[row*bw+k-i])));
 oh[gid]=z.x;ol[gid]=z.y;
}
'''


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def atomic(path: Path, value: object) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def dsconv(lib, a, b):
    """One batched DS convolution. a/b are (hi, lo, width) MPS tensors."""
    ah, al, aw = a
    bh, bl, bw = b
    n = ah.numel() // aw
    ow = aw + bw - 1
    oh = torch.empty(n * ow, dtype=torch.float32, device="mps")
    ol = torch.empty_like(oh)
    # MPS compile_shader obtains the grid from the first output tensor.
    lib.dsconv(oh, ol, ah, al, bh, bl, aw, bw)
    return oh, ol, ow


def split_poly(p: dict, width: int) -> tuple[np.ndarray, np.ndarray, int]:
    h = np.asarray(p["hi"], dtype=np.float32)
    l = np.asarray(p["lo"], dtype=np.float32)
    assert h.shape == l.shape and h.size <= width
    return np.pad(h, (0, width - h.size)), np.pad(l, (0, width - l.size)), h.size - 1


def decode(code: int, radices: list[int]) -> list[int]:
    ans = []
    for r in radices:
        ans.append(code % r)
        code //= r
    return ans


def selected_groups(blocks, E: int, digits: list[int]) -> tuple[int, list[dict]]:
    """Return the A-root bitmask and an auditable per-block selection record."""
    mask = 0
    records = []
    for block, d in zip(blocks, digits):
        option = block["options"][str(E)][d]
        if E < 12:
            picks = [block["lift_groups"][q] for q in option["selected_lift_positions"]]
            mode = "selected_lifts"
        elif option["mode"] == "original":
            picks, mode = [block["original_group"]], "original"
        else:
            picks, mode = list(block["lift_groups"]), "all_lifts"
        mask |= sum(1 << int(x) for x in picks)
        records.append({"block": int(block["block"]), "option_index": int(d), "mode": mode,
                        "A_group_indices": [int(x) for x in picks]})
    return mask, records


def roots_reconstruct(groups: list[np.ndarray], mask: int) -> tuple[np.ndarray, np.ndarray]:
    """Independent float64 product of individual root factors, not blocks."""
    a = np.array([1.0])
    b = np.array([1.0])
    for j, factor in enumerate(groups):
        if (mask >> j) & 1:
            a = np.convolve(a, factor)
        else:
            b = np.convolve(b, factor)
    a /= a[0]
    b /= b[0]
    # The two endpoints equal one by the already-proved root-balance identity.
    a[0] = a[-1] = 1.0
    b[0] = b[-1] = 1.0
    return a, b


def gpu_poly(h: np.ndarray, l: np.ndarray, row: int, width: int, degree: int) -> np.ndarray:
    v = h[row * width:(row + 1) * width].astype(np.float64) + l[row * width:(row + 1) * width].astype(np.float64)
    v = v[:degree + 1] / v[0]
    v[0] = v[-1] = 1.0
    return v


def choose_options(data: dict, E: int):
    blocks = data["blocks"]
    opts = [list(b["options"][str(E)]) for b in blocks]
    if E == 12:
        # Root-balance convention: first parent remains original; the other
        # two independently choose original/all-lifts.
        opts[0] = [opts[0][0]]
    return opts


def run_batch(lib, options, codes: np.ndarray, eta: float):
    n = codes.size
    radices = [len(x) for x in options]
    digits = np.empty((n, 3), dtype=np.int64)
    q = codes.copy()
    for j, r in enumerate(radices):
        digits[:, j] = q % r
        q //= r
    max_a = [max(len(x["Ablock"]["hi"]) for x in o) for o in options]
    max_b = [max(len(x["Bblock"]["hi"]) for x in o) for o in options]
    packed = []
    degrees_a = np.empty(n, dtype=np.uint32)
    degrees_b = np.empty(n, dtype=np.uint32)
    for side, widths in (("Ablock", max_a), ("Bblock", max_b)):
        hs, ls = [], []
        for j, (choices, width) in enumerate(zip(options, widths)):
            h = np.zeros((n, width), dtype=np.float32)
            l = np.zeros_like(h)
            for i in range(n):
                ph, pl, deg = split_poly(choices[int(digits[i, j])][side], width)
                h[i], l[i] = ph, pl
                if side == "Ablock":
                    degrees_a[i] = degrees_a[i] + deg if j else deg
                else:
                    degrees_b[i] = degrees_b[i] + deg if j else deg
            hs.append(torch.from_numpy(h.ravel()).to("mps"))
            ls.append(torch.from_numpy(l.ravel()).to("mps"))
        x01 = dsconv(lib, (hs[0], ls[0], widths[0]), (hs[1], ls[1], widths[1]))
        x = dsconv(lib, x01, (hs[2], ls[2], widths[2]))
        packed.append(x)
        del hs, ls, x01
    ah, al, aw = packed[0]
    bh, bl, bw = packed[1]
    # The four products are DS in Metal.  The final MPS vector screen adds
    # hi+lo and normalizes in float32; retained rows get a separate CPU64
    # reconstruction from individual root factors before they are reported.
    da = torch.from_numpy(degrees_a.astype(np.int64)).to("mps")
    db = torch.from_numpy(degrees_b.astype(np.int64)).to("mps")
    av = (ah + al).reshape(n, aw); bv = (bh + bl).reshape(n, bw)
    av = av / av[:, :1]; bv = bv / bv[:, :1]
    ia = torch.arange(aw, device="mps")[None, :]
    ib = torch.arange(bw, device="mps")[None, :]
    active_a = (ia > 0) & (ia < da[:, None])
    active_b = (ib > 0) & (ib < db[:, None])
    active = torch.cat((active_a, active_b), dim=1)
    values = torch.cat((av, bv), dim=1)
    score = (torch.relu(-values) + torch.relu(values - 1.0)).masked_fill(~active, 0.0).sum(dim=1)
    guarded = torch.where(active, values, torch.ones_like(values))
    mn, mx = guarded.min(dim=1).values, guarded.max(dim=1).values
    ma = torch.where(active, values.abs(), torch.ones_like(values)).max(dim=1).values
    frac = ((values > eta) & (values < 1.0 - eta) & active).any(dim=1)
    torch.mps.synchronize()
    result = {
        "digits": digits, "degree_A": degrees_a, "degree_B": degrees_b, "width_A": aw, "width_B": bw,
        "score": score.cpu().numpy(), "min": mn.cpu().numpy(), "max": mx.cpu().numpy(),
        "maxabs": ma.cpu().numpy(), "fractional": frac.cpu().numpy().astype(bool),
        "Ahi": ah.cpu().numpy(), "Alo": al.cpu().numpy(), "Bhi": bh.cpu().numpy(), "Blo": bl.cpu().numpy(),
    }
    del ah, al, bh, bl, da, db, av, bv, ia, ib, active_a, active_b, active, values, guarded, score, mn, mx, ma, frac, packed
    return result


def candidate_record(data, E: int, code: int, digits: list[int], gpu: dict, i: int) -> dict:
    mask, mapping = selected_groups(data["blocks"], E, digits)
    return {"code": int(code), "mixed_radix_digits_block0_first": [int(x) for x in digits],
            "A_root_group_bitmask_decimal": str(mask), "A_root_group_bitmask_hex": hex(mask),
            "mask_mapping": mapping, "degree_A": int(gpu["degree_A"][i]), "degree_B": int(gpu["degree_B"][i]),
            "gpu_ds_violation_l1": float(gpu["score"][i]), "gpu_ds_min_interior_or_endpoint": float(gpu["min"][i]),
            "gpu_ds_max_interior_or_endpoint": float(gpu["max"][i]), "gpu_ds_max_abs": float(gpu["maxabs"][i]),
            "gpu_ds_has_fractionality": bool(gpu["fractional"][i])}


def recheck_top(entry: dict, groups: list[np.ndarray], C: np.ndarray) -> None:
    for row in entry["top64"]:
        mask = int(row["A_root_group_bitmask_decimal"])
        a, b = roots_reconstruct(groups, mask)
        prod = np.convolve(a, b)
        target = C / C[0]
        row["cpu_float64_independent_root_reconstruction"] = {
            "A_coefficients_ascending": a.tolist(), "B_coefficients_ascending": b.tolist(),
            "convolution_max_abs_residual": float(np.max(np.abs(prod - target))),
            "minimum_actual_coefficient": float(min(a.min(), b.min())),
            "maximum_actual_coefficient": float(max(a.max(), b.max())),
            "has_fractionality_eta_1e-5": bool(np.any((a[1:-1] > 1e-5) & (a[1:-1] < 1 - 1e-5)) or np.any((b[1:-1] > 1e-5) & (b[1:-1] < 1 - 1e-5))),
        }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--options", type=Path, required=True)
    ap.add_argument("--roots-input", type=Path, required=True)
    ap.add_argument("--candidate-index", type=int, default=2)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=32768)
    ap.add_argument("--max-codes", type=int)
    ap.add_argument("--smoke", action="store_true", help="run E=1,6,12, each at <=1024 codes")
    ap.add_argument("--eta", type=float, default=1e-5)
    a = ap.parse_args()
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data is not mounted; refusing to write elsewhere")
    if a.batch_size <= 0 or a.batch_size > 131072:
        raise SystemExit("batch size must be 1..131072")
    if not torch.backends.mps.is_available():
        raise SystemExit("MPS unavailable")
    data = json.loads(a.options.read_text())
    raw = json.loads(a.roots_input.read_text())["candidates"][a.candidate_index]
    if data["candidate"] != raw["name"] or data["L"] != 12 or len(data["blocks"]) != 3:
        raise SystemExit("expected prepared G L=12 candidate-index 2 with three blocks")
    groups = [np.asarray(g["coefficients_ascending"], dtype=np.float64) for g in raw["groups"]]
    C = np.asarray(raw["C_coefficients_ascending"], dtype=np.float64)
    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    script = Path(__file__)
    report = {"format": "streaming-balanced-double-single-mps-v1", "status": "running",
              "source_sha256": sha(script), "options_sha256": sha(a.options), "roots_input_sha256": sha(a.roots_input),
              "candidate": data["candidate"], "candidate_index": a.candidate_index,
              "parameters": {"batch_size": a.batch_size, "eta": a.eta, "max_codes": a.max_codes, "smoke": a.smoke},
              "algorithm": {"device": "mps", "float_format": "double-single float32 hi+lo", "local_block_count": 3,
                 "batch_DS_convolutions": 4, "mixed_radix": "block 0 least significant", "E_L_rule": "block0 original; blocks1,2 original/all_lifts",
                 "endpoint_normalization": "true endpoints set to 1 using proved root-balance identity; padded tails excluded"},
              "versions": {"python": sys.version, "torch": torch.__version__, "platform": platform.platform()},
              "E_progress": [], "claim_boundary": "MPS DS ranks numerical candidates only; it is not a proof of acceptance or exclusion."}
    # A normal interrupted scan resumes its contiguous coverage.  Smoke files
    # deliberately do not resume into a full scan because their parameters
    # describe a different, bounded experiment.
    prior: dict[int, dict] = {}
    if a.checkpoint.exists():
        old = json.loads(a.checkpoint.read_text())
        compatible = all(old.get(k) == report.get(k) for k in ("format", "source_sha256", "options_sha256", "roots_input_sha256", "candidate", "candidate_index"))
        if compatible and old.get("parameters", {}).get("smoke") == a.smoke:
            prior = {int(x["E"]): x for x in old.get("E_progress", [])}
            report["resumed_from_checkpoint"] = str(a.checkpoint)
    started = perf_counter()
    try:
        lib = torch.mps.compile_shader(METAL)
        Es = [1, 6, 12] if a.smoke else list(range(1, 13))
        for E in Es:
            options = choose_options(data, E)
            rad = [len(x) for x in options]
            total = math.prod(rad)
            end = min(total, 1024 if a.smoke else (a.max_codes if a.max_codes is not None else total))
            old_entry = prior.get(E, {})
            old_coverage = old_entry.get("coverage_completed_half_open", [0, 0])
            start = int(old_coverage[1]) if old_coverage and old_coverage[0] == 0 else 0
            start = min(start, end)
            best = list(old_entry.get("top64", []))
            counts = dict(old_entry.get("counts", {"codes_seen": 0, "fractional_masks": 0, "closed_unit_interval_fractional_masks": 0}))
            smoke_errors = []
            fair = None
            entry = {"E": E, "status": "running", "radices_block0_first": rad, "total_codes": total,
                     "coverage_completed_half_open": [0, start], "counts": counts, "top64": best}
            for lo in range(start, end, a.batch_size):
                hi = min(end, lo + a.batch_size)
                codes = np.arange(lo, hi, dtype=np.int64)
                gpu = run_batch(lib, options, codes, a.eta)
                counts["fractional_masks"] += int(gpu["fractional"].sum())
                counts["closed_unit_interval_fractional_masks"] += int((gpu["fractional"] & (gpu["min"] >= 0) & (gpu["max"] <= 1)).sum())
                # Retain only the batch's best 64 before materialising costly
                # root masks/mappings.  This keeps resident memory O(batch).
                pick = np.lexsort((codes, gpu["score"]))[:min(64, len(codes))]
                best.extend(candidate_record(data, E, int(codes[i]), gpu["digits"][i].tolist(), gpu, int(i)) for i in pick)
                best = sorted(best, key=lambda z: (z["gpu_ds_violation_l1"], z["code"]))[:64]
                counts["codes_seen"] += len(codes)
                # Every smoke code is independently reconstructed from 39 root factors.
                if a.smoke:
                    for i, code in enumerate(codes):
                        mask, _ = selected_groups(data["blocks"], E, gpu["digits"][i].tolist())
                        ca, cb = roots_reconstruct(groups, mask)
                        ga = gpu_poly(gpu["Ahi"], gpu["Alo"], i, gpu["width_A"], int(gpu["degree_A"][i]))
                        gb = gpu_poly(gpu["Bhi"], gpu["Blo"], i, gpu["width_B"], int(gpu["degree_B"][i]))
                        smoke_errors.append(max(float(np.max(np.abs(ga - ca))), float(np.max(np.abs(gb - cb)))))
                        if E == 12 and int(code) == 0:
                            H = np.asarray(raw["H_coefficients_ascending"], dtype=np.float64)
                            HL = np.zeros((len(H) - 1) * 12 + 1); HL[::12] = H
                            fair = {"code": 0, "A_vs_exact_integer_H_max_abs": float(np.max(np.abs(ga - H))),
                                    "B_vs_exact_integer_H_of_x12_max_abs": float(np.max(np.abs(gb - HL))),
                                    "cpu_root_A_vs_H_max_abs": float(np.max(np.abs(ca - H))),
                                    "cpu_root_B_vs_H_of_x12_max_abs": float(np.max(np.abs(cb - HL)))}
                entry = {"E": E, "status": "running", "radices_block0_first": rad, "total_codes": total,
                         "coverage_completed_half_open": [0, hi], "counts": counts, "top64": best}
                report["E_progress"] = [z for z in report["E_progress"] if z["E"] != E] + [entry]
                report["elapsed_seconds"] = perf_counter() - started
                atomic(a.checkpoint, report)
                del gpu
                torch.mps.empty_cache()
            entry["status"] = "complete" if end == total else "partial_unknown"
            entry["coverage_completed_half_open"] = [0, end]
            if a.smoke:
                entry["smoke_gpu_vs_independent_cpu_float64_root_reconstruction"] = {"codes_checked": end,
                    "max_abs_error": float(max(smoke_errors, default=0.0)), "per_code_errors": smoke_errors}
                if fair is not None:
                    entry["fair_E12_exact_integer_coefficient_check"] = fair
            recheck_top(entry, groups, C)
            report["E_progress"] = [z for z in report["E_progress"] if z["E"] != E] + [entry]
            report["elapsed_seconds"] = perf_counter() - started
            atomic(a.checkpoint, report)
        report["status"] = "smoke_complete" if a.smoke else ("complete" if all(x["status"] == "complete" for x in report["E_progress"]) else "partial_unknown")
        report["elapsed_seconds"] = perf_counter() - started
        atomic(a.checkpoint, report)
        print(json.dumps({"status": report["status"], "checkpoint": str(a.checkpoint), "E": Es}, indent=2))
    finally:
        if torch.backends.mps.is_available():
            torch.mps.synchronize()
            torch.mps.empty_cache()


if __name__ == "__main__":
    main()
