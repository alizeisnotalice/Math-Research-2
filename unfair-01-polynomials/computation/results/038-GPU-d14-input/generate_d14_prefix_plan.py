#!/usr/bin/env python3
"""Generate the d=14 necessary-prefix pool and a deterministic start-row plan.

This program deliberately does not import the GPU worker or its prefix module.
All prefix arithmetic below is integer arithmetic.  It prepares input only;
it neither initializes MPS nor makes any mathematical exclusion claim.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

PAIR = ((0, 0), (1, 0), (0, 1))


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def decode(code: str) -> tuple[list[int], list[int]]:
    a, b = [1], [1]
    for char in code:
        ai, bi = PAIR[ord(char) - ord("0")]
        a.append(ai)
        b.append(bi)
    return a, b


def convolution(a: list[int], b: list[int], end: int) -> list[int]:
    return [sum(a[i] * b[k - i]
                for i in range(max(0, k - len(b) + 1), min(k, len(a) - 1) + 1))
            for k in range(end + 1)]


def legal(code: str, d: int) -> bool:
    a, b = decode(code)
    return (set(convolution(a, b, d - 1)) <= {0, 1}
            and sum(a[i] * b[d - i] for i in range(1, d)) == 0)


def reasons(code: str, d: int) -> list[str]:
    a, b = decode(code)
    product = convolution(a, b, 2 * d - 2)
    answer: list[str] = []
    if any(value > 1 for value in product):
        answer.append("prefix_product_coefficient_gt_1")
    for i in range(1, d):
        # The coefficient past degree 2d-2 is exactly zero, not an error.
        baseline = product[d + i] if d + i < len(product) else 0
        if baseline >= 1 and a[i] + b[i] > 0:
            answer.append("forced_positive_fractional_cross_at_d_plus_i")
            break
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--d", type=int, default=14)
    parser.add_argument("--minimum-starts", type=int, default=8192)
    parser.add_argument("--m", type=int, default=60)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--margin", type=float, default=.001)
    parser.add_argument("--seed", type=int, default=2026111615)
    args = parser.parse_args()
    assert args.d == 14 and args.minimum_starts >= 8192 and args.minimum_starts % 8192 == 0
    assert args.m > args.d and args.n > args.d and 0 < args.margin < .5

    raw = ["".join(chars) for chars in itertools.product("012", repeat=args.d - 1)]
    valid = [code for code in raw if legal(code, args.d)]
    survivors = [code for code in valid if not reasons(code, args.d)]
    assert survivors
    planned_starts = max(args.minimum_starts, ((len(survivors) + 8191) // 8192) * 8192)
    assignments = [survivors[row % len(survivors)] for row in range(planned_starts)]
    counts = Counter(assignments)
    assert set(counts) == set(survivors)
    assert sum(counts.values()) == planned_starts
    filter_counts = {reason: sum(reason in reasons(code, args.d) for code in valid)
                     for reason in ("prefix_product_coefficient_gt_1", "forced_positive_fractional_cross_at_d_plus_i")}
    output = {
        "schema": "d14-necessary-prefix-plan-v1",
        "status": "cpu_generated_not_gpu_started",
        "model": "A0=B0=1; i<d pairs in {(0,0),(1,0),(0,1)}; low convolution binary; d-cross zero; listed filters necessary only",
        "parameters": {"d": args.d, "m": args.m, "n": args.n, "margin": args.margin,
                       "seed": args.seed, "minimum_starts": args.minimum_starts,
                       "coverage_plan_starts": planned_starts},
        "raw_ternary_count": len(raw),
        "valid_low_prefix_count": len(valid),
        "survivor_count": len(survivors),
        "valid_low_prefix_codes_sha256": digest(valid),
        "survivor_codes_sha256": digest(survivors),
        "strict_filter_rejected_by_reason_counts": filter_counts,
        "survivor_codes": survivors,
        "coverage_assignment_codes": assignments,
        "per_prefix_assignment_counts": [{"code": code, "count": counts[code]} for code in survivors],
        "minimum_starts_covers_pool": args.minimum_starts >= len(survivors),
        "coverage_plan_is_8192_multiple": planned_starts % 8192 == 0,
        "coverage_plan_covers_every_survivor": True,
        "claim_boundary": "Exact finite prefix generation only. Neither filter survival nor any later numerical MPS trajectory is an exclusion theorem.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: output[key] for key in ("raw_ternary_count", "valid_low_prefix_count", "survivor_count", "minimum_starts_covers_pool")}, sort_keys=True))


if __name__ == "__main__":
    main()
