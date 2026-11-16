#!/usr/bin/env python3
"""Bounded exact Q[x] factor-allocation search for P31 at (p,q)=(3,5).

The search domain is deliberately finite: binary, monic, constant-one H of
degrees 18..24 whose H(x^3)H(x^5) is binary.  It uses a support DFS and
exact integer arithmetic only.  No conclusion about arbitrary real factors
or degrees outside this interval is licensed by an empty result.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
import sys
import time
from collections import Counter
from pathlib import Path

import sympy as sp

P, Q = 3, 5
X = sp.symbols("x")


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_json(value: object) -> str:
    return sha_bytes(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def atom(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def conv(a: list[int], b: list[int]) -> list[int]:
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj
    return result


def stretch(a: list[int], scale: int) -> list[int]:
    result = [0] * ((len(a) - 1) * scale + 1)
    for i, ai in enumerate(a):
        result[scale * i] = ai
    return result


def first_negative(a: list[int]) -> dict[str, object] | None:
    for degree, coefficient in enumerate(a):
        if coefficient < 0:
            return {"degree": degree, "coefficient": str(coefficient)}
    return None


def nonnegative(a: list[int]) -> bool:
    return all(coefficient >= 0 for coefficient in a)


def binary_stretch_product_is_binary(h: list[int]) -> bool:
    c = conv(stretch(h, P), stretch(h, Q))
    return all(coefficient in (0, 1) for coefficient in c)


def support_dfs(max_degree: int):
    """Yield all supports S with 0 in S and binary 3S+5S sumset.

    Appending t adds precisely {3t+5s,3s+5t:s in S} union {8t}; it is
    accepted only if these 2|S|+1 degrees are mutually new and avoid the old
    sumset.  Thus every yielded support has a binary H(x^3)H(x^5).
    """
    def descend(support: tuple[int, ...], sumset: frozenset[int]):
        yield support
        for t in range(support[-1] + 1, max_degree + 1):
            new = {3 * t + 5 * s for s in support}
            new.update(3 * s + 5 * t for s in support)
            new.add(8 * t)
            if len(new) != 2 * len(support) + 1 or not sumset.isdisjoint(new):
                continue
            yield from descend(support + (t,), sumset | frozenset(new))

    yield from descend((0,), frozenset({0}))


def h_from_support(support: tuple[int, ...]) -> list[int]:
    h = [0] * (support[-1] + 1)
    for exponent in support:
        h[exponent] = 1
    return h


def coefficient_list(poly: sp.Poly) -> list[int]:
    values = [poly.nth(i) for i in range(poly.degree() + 1)]
    if any(value.q != 1 for value in values):
        raise AssertionError("monic rational factor did not have integer coefficients")
    return [int(value) for value in values]


def factor_h(h: list[int]) -> tuple[list[sp.Poly], list[int]]:
    poly = sp.Poly.from_list(list(reversed(h)), X, domain=sp.ZZ)
    content, factor_data = sp.factor_list(poly.as_expr(), X)
    if content != 1:
        raise AssertionError("unexpected nonunit H content")
    factors = [sp.Poly(factor, X, domain=sp.QQ) for factor, _ in factor_data]
    multiplicities = [int(mult) for _, mult in factor_data]
    if not all(factor.LC() == 1 and factor.is_irreducible for factor in factors):
        raise AssertionError("factor_list did not return monic irreducible QQ factors")
    rebuilt = [1]
    for factor, multiplicity in zip(factors, multiplicities):
        factor_coefficients = coefficient_list(factor)
        for _ in range(multiplicity):
            rebuilt = conv(rebuilt, factor_coefficients)
    if rebuilt != h:
        raise AssertionError("exact factor product mismatch")
    return factors, multiplicities


def multiply_factor_powers(factors: list[list[int]], exponents: tuple[int, ...]) -> list[int]:
    result = [1]
    for factor, exponent in zip(factors, exponents):
        for _ in range(exponent):
            result = conv(result, factor)
    return result


def append_jsonl(path: Path, value: object) -> None:
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def checkpoint(out: Path, last_index: int, counts: Counter, started: float) -> None:
    atom(out / "checkpoint.json", {
        "status": "running",
        "last_completed_h_index": last_index,
        "counts": dict(counts),
        "elapsed_seconds": time.monotonic() - started,
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-degree", type=int, default=18)
    parser.add_argument("--max-degree", type=int, default=24)
    parser.add_argument("--checkpoint-every", type=int, default=50)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not (18 <= args.min_degree <= args.max_degree <= 24):
        raise SystemExit("this bounded run is defined only for degrees 18..24")
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data is not mounted")
    if args.out.exists() and any(args.out.iterdir()) and not args.resume:
        raise SystemExit("refusing nonempty output without --resume")
    args.out.mkdir(parents=True, exist_ok=True)
    coverage_path, allocation_path = args.out / "h-coverage.jsonl", args.out / "signed-allocations.jsonl"
    if args.resume:
        saved = json.loads((args.out / "checkpoint.json").read_text())
        start_index = int(saved["last_completed_h_index"]) + 1
        counts = Counter(saved["counts"])
    else:
        start_index, counts = 0, Counter()
        for path in (coverage_path, allocation_path):
            path.touch(exist_ok=False)

    all_supports = [support for support in support_dfs(args.max_degree)
                    if args.min_degree <= support[-1] <= args.max_degree]
    degree_counts = Counter(support[-1] for support in all_supports)
    started = time.monotonic()
    metadata = {
        "scope": {"p": P, "q": Q, "degrees_inclusive": [args.min_degree, args.max_degree],
                  "H": "all support-DFS binary monic constant-one H with binary H(x^3)H(x^5)",
                  "factors": "all allocations of SymPy exact QQ irreducible factor multiplicities"},
        "enumerator": {"rule": "append t if new crossed degrees have size 2|S|+1 and avoid old 3S+5S",
                       "actual_count_by_degree": {str(d): degree_counts[d] for d in range(args.min_degree, args.max_degree + 1)},
                       "total_H": len(all_supports)},
        "exact_arithmetic": "Python integer convolutions; all QQ monic factors are asserted integral before allocation products",
        "factorization_assumption": "SymPy factor_list over QQ and Poly.is_irreducible are treated as the exact algebra backend; this run does not supply an independent formal irreducibility certificate.",
        "software": {"python": sys.version, "sympy": sp.__version__, "platform": platform.platform()},
        "source_sha256": sha_bytes(Path(__file__).read_bytes()),
    }
    atom(args.out / "metadata.json", metadata)

    for index, support in enumerate(all_supports):
        if index < start_index:
            continue
        h = h_from_support(support)
        if not binary_stretch_product_is_binary(h):
            raise AssertionError("DFS emitted a nonbinary C")
        factors, multiplicities = factor_h(h)
        factor_coefficients = [coefficient_list(factor) for factor in factors]
        allocation_counts = Counter()
        signed_ids: list[str] = []
        h_id = f"{index:05d}-d{support[-1]}"
        for exponents in itertools.product(*(range(multiplicity + 1) for multiplicity in multiplicities)):
            allocation_counts["all"] += 1
            f = multiply_factor_powers(factor_coefficients, exponents)
            complement = tuple(multiplicity - exponent for multiplicity, exponent in zip(multiplicities, exponents))
            g = multiply_factor_powers(factor_coefficients, complement)
            if conv(f, g) != h:
                raise AssertionError("allocation factor product mismatch")
            if f[0] <= 0 or g[0] <= 0:
                allocation_counts["nonpositive_constant"] += 1
                continue
            allocation_counts["positive_constant"] += 1
            if nonnegative(f) and nonnegative(g):
                allocation_counts["both_factor_nonnegative_skipped"] += 1
                continue
            allocation_counts["signed_checked"] += 1
            a = conv(stretch(f, P), stretch(g, Q))
            b = conv(stretch(g, P), stretch(f, Q))
            first_a, first_b = first_negative(a), first_negative(b)
            allocation_id = h_id + ":" + ",".join(map(str, exponents))
            signed_ids.append(allocation_id)
            record = {
                "allocation_id": allocation_id, "h_id": h_id, "degree": support[-1], "support": list(support),
                "H": h, "factor_coefficients_low_to_high": factor_coefficients, "multiplicities": multiplicities,
                "allocation_exponents_for_F": list(exponents), "F": f, "G": g,
                "first_negative_A_equals_F_x3_G_x5": first_a,
                "first_negative_B_equals_G_x3_F_x5": first_b,
            }
            if first_a is None and first_b is None:
                allocation_counts["both_cross_nonnegative"] += 1
                record["status"] = "candidate_both_cross_nonnegative"
                record["A"] = a
                record["B"] = b
            else:
                allocation_counts["cross_negative_witnessed"] += 1
                record["status"] = "rejected_by_exact_cross_negative_witness"
            append_jsonl(allocation_path, record)
        coverage = {
            "h_id": h_id, "degree": support[-1], "support": list(support), "H": h,
            "h_sha256": sha_json(h), "factor_coefficients_low_to_high": factor_coefficients,
            "multiplicities": multiplicities, "factor_count": len(factors),
            "factorization_verified_by_exact_integer_product": True,
            "allocation_counts": dict(allocation_counts), "signed_allocation_ids": signed_ids,
        }
        append_jsonl(coverage_path, coverage)
        counts.update({"H_completed": 1})
        counts.update({"allocations_all": allocation_counts["all"],
                       "allocations_positive_constant": allocation_counts["positive_constant"],
                       "allocations_nonpositive_constant": allocation_counts["nonpositive_constant"],
                       "allocations_nonnegative_skipped": allocation_counts["both_factor_nonnegative_skipped"],
                       "signed_allocations_checked": allocation_counts["signed_checked"],
                       "cross_negative_witnessed": allocation_counts["cross_negative_witnessed"],
                       "both_cross_nonnegative_candidates": allocation_counts["both_cross_nonnegative"]})
        if (index + 1) % args.checkpoint_every == 0:
            checkpoint(args.out, index, counts, started)

    final = {**metadata, "status": "complete", "elapsed_seconds": time.monotonic() - started,
             "counts": dict(counts), "output_hashes": {name: sha_bytes((args.out / name).read_bytes())
             for name in ("h-coverage.jsonl", "signed-allocations.jsonl")},
             "claim_boundary": "This is a bounded exact rational factor-allocation search only. It does not cover arbitrary real factors, degrees outside 18..24, or prove the proposed positivity descent statement."}
    atom(args.out / "final.json", final)
    atom(args.out / "checkpoint.json", {"status": "complete", "last_completed_h_index": len(all_supports) - 1,
                                         "counts": dict(counts), "elapsed_seconds": final["elapsed_seconds"]})
    print(json.dumps({"status": final["status"], "counts": final["counts"], "final_sha256": sha_bytes((args.out / "final.json").read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
