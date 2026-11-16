#!/usr/bin/env python3
"""Independent integer audit for a d=13 necessary-prefix plan.

This verifier intentionally neither imports nor executes the generator or the
GPU worker. It recomputes the full finite pool and checks every planned row.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

PAIR = {"0": (0, 0), "1": (1, 0), "2": (0, 1)}


def canonical_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def expand(code: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    aa, bb = [1], [1]
    for symbol in code:
        ai, bi = PAIR[symbol]
        aa.append(ai)
        bb.append(bi)
    return tuple(aa), tuple(bb)


def coeff(aa: tuple[int, ...], bb: tuple[int, ...], k: int) -> int:
    return sum(aa[i] * bb[k - i]
               for i in range(max(0, k - len(bb) + 1), min(k, len(aa) - 1) + 1))


def allowed(code: str, d: int) -> bool:
    aa, bb = expand(code)
    return (all(coeff(aa, bb, k) in (0, 1) for k in range(d))
            and sum(aa[i] * bb[d - i] for i in range(1, d)) == 0)


def reject_reasons(code: str, d: int) -> tuple[str, ...]:
    aa, bb = expand(code)
    product = [coeff(aa, bb, k) for k in range(2 * d - 1)]
    results: list[str] = []
    if max(product) > 1:
        results.append("prefix_product_coefficient_gt_1")
    if any((product[d + i] if d + i < len(product) else 0) >= 1
           and aa[i] + bb[i] > 0 for i in range(1, d)):
        results.append("forced_positive_fractional_cross_at_d_plus_i")
    return tuple(results)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--prefix-module", type=Path, required=True)
    parser.add_argument("--worker-cpu-audit", type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.input.read_text())
    assert plan["schema"] == "d13-necessary-prefix-plan-v1"
    p = plan["parameters"]
    assert p == {"coverage_plan_starts": p["coverage_plan_starts"], "d": 13, "m": 56,
                 "margin": .001, "n": 96, "requested_starts": 32768, "seed": 2026111611}
    d = p["d"]
    raw = ["".join(chars) for chars in itertools.product("012", repeat=d - 1)]
    valid = [code for code in raw if allowed(code, d)]
    survivors = [code for code in valid if not reject_reasons(code, d)]
    listed = plan["survivor_codes"]
    assert len(listed) == len(set(listed)), "duplicate survivor code"
    assert listed == survivors, "generator order or membership differs from independent enumeration"
    assert plan["raw_ternary_count"] == len(raw) == 3 ** 12
    assert plan["valid_low_prefix_count"] == len(valid)
    assert plan["survivor_count"] == len(survivors)
    assert plan["valid_low_prefix_codes_sha256"] == canonical_hash(valid)
    assert plan["survivor_codes_sha256"] == canonical_hash(survivors)
    for code in survivors:
        assert allowed(code, d) and not reject_reasons(code, d)
    assignments = plan["coverage_assignment_codes"]
    assert len(assignments) == p["coverage_plan_starts"]
    assert assignments == [survivors[index % len(survivors)] for index in range(len(assignments))]
    counts = Counter(assignments)
    listed_counts = plan["per_prefix_assignment_counts"]
    assert [row["code"] for row in listed_counts] == survivors
    assert [row["count"] for row in listed_counts] == [counts[code] for code in survivors]
    assert set(counts) == set(survivors)
    assert sum(counts.values()) == len(assignments)
    q, r = divmod(len(assignments), len(survivors))
    expected_counts = [q + (index < r) for index in range(len(survivors))]
    assert [counts[code] for code in survivors] == expected_counts
    assert plan["requested_starts_covers_pool"] == (32768 >= len(survivors))
    assert plan["coverage_plan_covers_every_survivor"] is True
    assert p["coverage_plan_starts"] == max(32768, len(survivors))
    assert "1..13" in args.worker.read_text() and "d <= 13" in args.worker.read_text()
    worker_audit = json.loads(args.worker_cpu_audit.read_text())
    worker_fields = {
        "d": worker_audit["d"],
        "raw_ternary_count": worker_audit["raw_ternary_code_count"],
        "valid_low_prefix_count": worker_audit["valid_low_prefix_count"],
        "survivor_count": worker_audit["survivor_count"],
        "valid_low_prefix_codes_sha256": worker_audit["valid_low_prefix_codes_sha256"],
        "survivor_codes_sha256": worker_audit["survivor_codes_sha256"],
        "filter_rejections": worker_audit["strict_filter_rejected_by_reason_counts"],
    }
    independently_recomputed = {
        "d": d,
        "raw_ternary_count": len(raw),
        "valid_low_prefix_count": len(valid),
        "survivor_count": len(survivors),
        "valid_low_prefix_codes_sha256": canonical_hash(valid),
        "survivor_codes_sha256": canonical_hash(survivors),
        "filter_rejections": {
            reason: sum(reason in reject_reasons(code, d) for code in valid)
            for reason in ("prefix_product_coefficient_gt_1", "forced_positive_fractional_cross_at_d_plus_i")
        },
    }
    assert worker_fields == independently_recomputed, "worker CPU audit differs from independent arithmetic"
    assert worker_audit["all_survivors_cpu_audited"] is True
    assert worker_audit["float_arithmetic_used"] is False
    assert worker_audit["starts"] == 32768 and worker_audit["coverage_by_start_rows"] is True
    result = {
        "status": "passed_independent_integer_d13_plan_audit",
        "d": 13,
        "raw_ternary_count": len(raw),
        "valid_low_prefix_count": len(valid),
        "survivor_count": len(survivors),
        "valid_low_prefix_codes_sha256": canonical_hash(valid),
        "survivor_codes_sha256": canonical_hash(survivors),
        "strict_filter_rejected_by_reason_counts": independently_recomputed["filter_rejections"],
        "all_survivors_independently_checked": True,
        "all_coverage_rows_independently_checked": True,
        "per_prefix_assignment_count_min": min(counts.values()),
        "per_prefix_assignment_count_max": max(counts.values()),
        "coverage_plan_starts": len(assignments),
        "requested_starts": 32768,
        "requested_starts_covers_pool": 32768 >= len(survivors),
        "source_sha256": {"plan": hashlib.sha256(args.input.read_bytes()).hexdigest(),
                          "worker": hashlib.sha256(args.worker.read_bytes()).hexdigest(),
                          "prefix_module": hashlib.sha256(args.prefix_module.read_bytes()).hexdigest(),
                          "worker_cpu_audit": hashlib.sha256(args.worker_cpu_audit.read_bytes()).hexdigest(),
                          "verifier": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "claim_boundary": "Independent finite integer audit only; no MPS/GPU was initialized and it makes no factorization or exclusion claim.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
