#!/usr/bin/env python3
"""Check all first parity cuts of the twenty R11 sharp flat spectra."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time


CONSTANTS = tuple(map(Fraction, ("3/2", "2", "5/2", "3")))


def fwht(values):
    size = 1
    while size < len(values):
        for base in range(0, len(values), 2 * size):
            for offset in range(size):
                i, j = base + offset, base + offset + size
                first, second = values[i], values[j]
                values[i], values[j] = first + second, first - second
        size *= 2
    return values


def rank(points):
    basis = []
    for point in points:
        for vector in basis:
            point = min(point, point ^ vector)
        if point:
            basis.append(point)
            basis.sort(reverse=True)
    return len(basis)


def examine(row, index):
    support = frozenset(row["support"])
    coefficients = dict(zip(row["support"], row["coefficient_numerators_over_8"]))
    assert set(coefficients) == support and len(support) == 64
    assert set(map(abs, coefficients.values())) == {1}
    values = [coefficients.get(point, 0) for point in range(1 << 15)]
    outputs = Counter(fwht(values))
    assert set(outputs) <= {-8, 8}
    assert rank(support) == 15
    flag = 1 << 15
    assert rank(point | flag for point in support) == 15

    hist = Counter()
    successful = {str(constant): [] for constant in CONSTANTS}
    maximum_slack = None
    maximum_slack_cuts = []
    queries = {a ^ b for a, b in combinations(support, 2)}
    period_q = []
    for query in sorted(queries):
        pivot = 1 << (query.bit_length() - 1)
        reps = {point ^ (query if point & pivot else 0) for point in support}
        one = set()
        zero = set()
        for point in reps:
            left = coefficients.get(point, 0)
            right = coefficients.get(point ^ query, 0)
            if left + right:
                zero.add(point)
            if left - right:
                one.add(point)
        double = 64 - len(reps)
        assert len(one) == len(zero) == 64 - 3 * double // 2
        overlap = len(one & zero)
        assert overlap == 64 - 2 * double
        span_zero = rank(point | flag for point in zero)
        span_one = rank(point | flag for point in one)
        span_union = rank(point | flag for point in reps)
        assert span_union == 14
        intersection = span_zero + span_one - span_union
        if double == 32:
            period_q.append(query)
        hist[(double, len(one), overlap, intersection)] += 1
        slack = (intersection + 16) ** 2 - 16 * len(one)
        if maximum_slack is None or slack > maximum_slack:
            maximum_slack = slack
            maximum_slack_cuts = [query]
        elif slack == maximum_slack:
            maximum_slack_cuts.append(query)
        for constant in CONSTANTS:
            if (intersection + 8 * constant) ** 2 >= 4 * constant ** 2 * len(one):
                successful[str(constant)].append(query)

    assert sorted(row["periods"]) == [0] + period_q
    assert successful["2"] == period_q
    result = {
        "certificate_index": index,
        "original_row_index": row["row_index"],
        "original_source_sha256": row["source_sha256"],
        "nonreducing_cuts": (1 << 15) - 1 - len(queries),
        "strict_cuts": len(queries),
        "boolean_output_count": sum(outputs.values()),
        "periods": period_q,
        "successful_cut_counts": {key: len(value) for key, value in successful.items()},
        "successful_C2_masks": successful["2"],
        "max_C2_squared_slack": maximum_slack,
        "number_of_cuts_maximizing_slack": len(maximum_slack_cuts),
        "smallest_cut_maximizing_slack": min(maximum_slack_cuts),
        "rank_profiles": {",".join(map(str, key)): count for key, count in sorted(hist.items())},
    }
    if index == 0:
        result["witness_support"] = row["support"]
        result["witness_numerators_over_8"] = row["coefficient_numerators_over_8"]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    reference_bytes = args.reference.read_bytes()
    reference = json.loads(reference_bytes)
    rows = [examine(row, index) for index, row in enumerate(reference["rows"])]
    assert len(rows) == 20 and sum(not row["periods"] for row in rows) == 10
    result = {
        "schema": "a28-r048-sharp-aperiodic-cuts-v1",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference_sha256": sha256(reference_bytes).hexdigest(),
        "profile_order": ["occupied_double_cosets", "common_child_sparsity", "M_equal_child", "h"],
        "records": rows,
        "strict_cut_count": sum(row["strict_cuts"] for row in rows),
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers and fractions; no randomness",
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"verified all {result['strict_cut_count']} strict cuts in {len(rows)} sharp parents")


if __name__ == "__main__":
    main()
