#!/usr/bin/env python3
"""Exhaustive small Boolean spectra and equality structure for directional Fourier mass."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
import hashlib
import json
from pathlib import Path
import platform
import runpy
import time


def setup(helper):
    global fwht, indexes, parity, independent, solve, span, rref
    tools = runpy.run_path(str(helper))
    fwht, indexes, parity, independent, solve, span, rref = [
        tools[name] for name in ("fwht", "indexes", "parity", "independent", "solve", "span", "rref")
    ]


def direction_spec(n, rows):
    r = len(rows)
    R, M = 1 << r, 1 << n
    labs = tuple(sum(parity(q & u) << i for i, u in enumerate(rows)) for q in range(M))
    rest = []
    for i in range(n):
        q = 1 << i
        if len(independent(rows + rest + [q])) > r + len(rest):
            rest.append(q)
    assert len(rest) + r == n
    cosets = span(rest)
    D = span(rows)
    assert len({x ^ u for x in cosets for u in D}) == M
    return {"n": n, "r": r, "M": M, "R": R, "rows": rows, "labels": labs, "D": D, "cosets": cosets}


@cache
def label_rank(mask):
    return len(independent(indexes(mask)))


def equality_structure(config, truth, spectrum, band_masses):
    n, r, M, R = (config[k] for k in ("n", "r", "M", "R"))
    assert band_masses[1:] == [2 * (M // R)] * (R - 1)
    by_beta = [[] for _ in range(R)]
    for a, val in enumerate(spectrum):
        k = config["labels"][a]
        if val and k:
            by_beta[k].append(a)
    magnitudes = {abs(spectrum[a]) for points in by_beta for a in points}
    assert len(magnitudes) == 1
    a = next(iter(magnitudes))
    common = {x ^ by_beta[1][0] for x in by_beta[1]}
    dim = len(independent(common))
    assert len(common) == 1 << dim
    assert all(config["labels"][q] == 0 for q in common)
    assert all({x ^ points[0] for x in points} == common for points in by_beta[1:])
    alpha0 = by_beta[1][0] ^ by_beta[2][0] ^ by_beta[3][0]
    missing = {alpha0 ^ q for q in common}
    assert all(config["labels"][q] == 0 for q in missing)
    offzero = sorted(q for section in by_beta for q in section)
    whole = set(offzero) | missing
    assert len(whole) == 1 << (r + dim)
    assert len(independent(q ^ min(whole) for q in whole)) == r + dim

    phase = solve([q | M for q in offzero], [int(spectrum[q] < 0) for q in offzero])
    x0 = phase & (M - 1)
    signed_sum = sum(spectrum[q] * (1 - 2 * parity(q & x0)) for q in missing)
    assert abs(signed_sum) == (R // 2 - 1) * a * len(common)
    assert a * len(common) == 2 * (M // R)

    active = []
    majority = []
    minority = []
    for j, x in enumerate(config["cosets"]):
        row = [truth[x ^ u] for u in config["D"]]
        if len(set(row)) == 1:
            continue
        sign = 1 if sum(row) > 0 else -1
        minor = [i for i, s in enumerate(row) if s != sign]
        assert len(minor) == 1
        active.append(j)
        majority.append(int(sign < 0))
        minority.append(minor[0])
    assert len(active) * len(common) == M // R
    assert len(active) == 1 << len(independent(q ^ active[0] for q in active))
    augmented_quotient = [u | (M // R) for u in active]
    solve(augmented_quotient, majority)
    for b in range(r):
        solve(augmented_quotient, [(u >> b) & 1 for u in minority])
    assert any(x0 in {config["cosets"][j] ^ u for u in config["D"]} for j in active)
    return {
        "common_K_basis": independent(common), "common_K_order": len(common),
        "nonzero_labels": R - 1, "offzero_support": offzero, "missing_coset": sorted(missing),
        "Fourier_coefficient_abs_value": str(Fraction(a, M)),
        "Fourier_sign_character_physical_point": x0,
        "signed_target_average_abs": str(Fraction(abs(signed_sum), M * len(common))),
        "unsigned_target_average": str(Fraction(sum(abs(spectrum[q]) for q in missing), M * len(common))),
        "active_physical_cosets": active, "majority_negative_bits_on_active": majority,
        "minority_positions_on_active": minority,
    }


def check_case(cfg, truth, W, stats, witness):
    R, M, r = cfg["R"], cfg["M"], cfg["r"]
    masses = [0] * R
    for beta, val in zip(cfg["labels"], W):
        masses[beta] += abs(val)
    mask = sum((masses[b] != 0) << b for b in range(1, R))
    stats["function_direction_pairs"] += 1
    if label_rank(mask) < r:
        return
    stats["full_span_pairs"] += 1
    all_labels = mask.bit_count() == R - 1
    stats["full_span_with_missing_labels"] += not all_labels
    extra = sum(masses[1:]) - (2 * M - 2 * (M // R))
    assert extra >= 0
    if extra:
        assert 32 * extra >= 5 * M
        gap_units32 = (32 * extra, M)
        old = witness.get("minimum_positive_gap")
        if old is None or gap_units32[0] * old["scaled_gap_denominator"] < gap_units32[1] * old["scaled_gap_numerator"]:
            witness["minimum_positive_gap"] = {
                "scaled_gap_numerator": gap_units32[0], "scaled_gap_denominator": gap_units32[1],
                "gap": str(Fraction(extra, M)), "direction_basis": cfg["rows"],
                "truth_mask": hex(sum((v < 0) << j for j, v in enumerate(truth))),
                "nonzero_labels": mask.bit_count(),
            }
    else:
        exact = equality_structure(cfg, truth, W, masses)
        assert all_labels
        stats["equality_cases"] += 1
        stats["equality_with_nonflat_full_spectrum"] += len({abs(x) for x in W if x}) > 1
        stats["equality_where_unsigned_zero_band_average_fails"] += (
            exact["signed_target_average_abs"] != exact["unsigned_target_average"]
        )
        stats[f"equality_K_order_{exact['common_K_order']}"] += 1
        witness.setdefault("first_equality", {
            "direction_basis": cfg["rows"], "truth_mask": hex(sum((v < 0) << j for j, v in enumerate(truth))),
            "structure": exact,
        })


def job(shard):
    n, lo, hi = shard
    configs = [direction_spec(n, basis) for r in range(2, n + 1) for basis in rref(n, r)]
    stats = {r: Counter() for r in range(2, n + 1)}
    witnesses = {r: {} for r in stats}
    M = 1 << n
    for mask in range(lo, hi):
        truth = [1 - 2 * ((mask >> x) & 1) for x in range(M)]
        W = fwht(truth)
        for cfg in configs:
            r = cfg["r"]
            check_case(cfg, truth, W, stats[r], witnesses[r])
    return n, stats, witnesses


def explicit(n, r, truth):
    M = 1 << n
    spec = direction_spec(n, [1 << j for j in range(r)])
    W = fwht(truth)
    counts, records = Counter(), {}
    check_case(spec, truth, W, counts, records)
    A = sum(abs(val) for freq, val in enumerate(W) if freq & ((1 << r) - 1))
    return {
        "n": n, "r": r, "A": str(Fraction(A, M)),
        "Walsh_support_size": sum(val != 0 for val in W),
        "truth_mask": hex(sum((v < 0) << j for j, v in enumerate(truth))),
        "Walsh_histogram": dict(Counter(W)), "stats": dict(counts), "witness": records,
    }


def constructions():
    results = {}
    for r in (2, 3, 4):
        n = r + 3
        F = []
        for x in range(1 << n):
            u, w = x % (1 << r), x >> r
            in_active = (w & 1) and not parity(w & 6)
            minority = (1 << (r - 1)) ^ ((w >> 1) | ((w & 1) << (r - 1)))
            if in_active:
                sign = parity(w & 5) ^ int(u == minority)
            else:
                sign = ((w >> 1) & (w >> 2) & 1) ^ parity(w)
            F.append(1 - 2 * sign)
        out = explicit(n, r, F)
        assert out["stats"]["equality_cases"] == 1 and out["stats"]["equality_K_order_4"] == 1
        results[f"codimension_two_active_affine_rows_r_{r}"] = out
    F = [(-1) ** (((x >> 2) & 1) ^ ((x & 1) & ((x >> 1) & 1))) for x in range(8)]
    gap = explicit(3, 3, F)
    assert gap["A"] == "2" and gap["witness"]["minimum_positive_gap"]["gap"] == "1/4"
    assert gap["stats"]["full_span_with_missing_labels"] == 1
    results["twisted_quadratic_gap_one_quarter"] = gap
    negative = explicit(3, 3, [1 - 2 * (x & 1) for x in range(8)])
    assert not negative["stats"].get("full_span_pairs") and negative["A"] == "1"
    results["rank_deficient_character_boundary"] = negative
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("half_windows_helper", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    start = time.perf_counter()
    setup(args.half_windows_helper)
    tasks = [(2, 0, 16), (3, 0, 256)] + [(4, i, i + 2048) for i in range(0, 65536, 2048)]
    all_counts, all_witnesses = {}, {}
    with ProcessPoolExecutor(max_workers=args.workers, initializer=setup, initargs=(args.half_windows_helper,)) as pool:
        for n, counts, records in pool.map(job, tasks):
            for r, vals in counts.items():
                all_counts.setdefault((n, r), Counter()).update(vals)
                out = all_witnesses.setdefault((n, r), {})
                out.setdefault("first_equality", records[r].get("first_equality"))
                if out["first_equality"] is None and "first_equality" in records[r]:
                    out["first_equality"] = records[r]["first_equality"]
                candidate = records[r].get("minimum_positive_gap")
                if candidate is not None:
                    old = out.get("minimum_positive_gap")
                    if old is None or Fraction(candidate["gap"]) < Fraction(old["gap"]):
                        out["minimum_positive_gap"] = candidate
    result = {
        "python": platform.python_version(), "workers": args.workers,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "helper_source_sha256": hashlib.sha256(args.half_windows_helper.read_bytes()).hexdigest(),
        "status": "All rank-full cases satisfy the lower bound, structured equality, and positive-gap test.",
        "exhaustive": [{"n": n, "r": r, "counts": counts, "witnesses": all_witnesses[n, r]}
                       for (n, r), counts in sorted(all_counts.items())],
        "explicit_cases": constructions(), "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
