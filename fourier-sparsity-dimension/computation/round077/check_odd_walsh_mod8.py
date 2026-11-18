#!/usr/bin/env python3
"""Integer checks for the affine normalization and quadratic mod-8 residue."""

import argparse
from collections import Counter
from itertools import combinations
import json
import platform
import random
import time


SEED = 770828
EXHAUSTIVE_DIMENSIONS = (3, 4)
SAMPLE_COUNTS = {5: 2_048, 6: 1_024, 7: 512, 8: 256, 9: 128, 10: 64}


def parity_table(size):
    parity = [0] * size
    for x in range(1, size):
        parity[x] = parity[x >> 1] ^ (x & 1)
    return parity


def walsh_half(size, truth):
    values = [1 - 2 * ((truth >> x) & 1) for x in range(size)]
    stride = 1
    while stride < size:
        for start in range(0, size, 2 * stride):
            for x in range(start, start + stride):
                left, right = values[x], values[x + stride]
                values[x], values[x + stride] = left + right, left - right
        stride *= 2
    if any(value % 2 for value in values):
        raise AssertionError("Walsh coefficient not divisible by 2")
    return [value // 2 for value in values]


def algebraic_degree(values, axes, size):
    coefficients = values.copy()
    for axis in axes:
        for u in range(size):
            if u & axis:
                coefficients[u] ^= coefficients[u ^ axis]
    return max((bin(u).count("1") for u in range(size) if coefficients[u]), default=-1)


def check_table(n, truth, parity, axes, pairs, direct, audit, digest):
    size = 1 << n
    support = [x for x in range(size) if truth >> x & 1]
    weight = len(support)
    xor_support = 0
    for x in support:
        xor_support ^= x
    global_sign = 1 - 2 * (((weight + 1) // 2) % 2)
    centered_support = [x ^ xor_support for x in support]
    first_moments = [sum(bool(x & axis) for x in centered_support) for axis in axes]
    if any(moment % 2 for moment in first_moments):
        raise AssertionError((n, truth, "odd centered first moment"))
    linear_axes = [axis for axis, moment in zip(axes, first_moments) if moment % 4 == 2]
    quadratic_pairs = [pair for pair in pairs if sum(x & pair == pair for x in centered_support) % 2]
    constant = ((-global_sign * weight - 1) // 4) % 2
    linear_mask = sum(linear_axes)
    quadratic_pattern = 0
    raw_values, corrected_values = [], []
    for u, half in enumerate(walsh_half(size, truth)):
        zero_term = size // 2 if u == 0 else 0
        if direct:
            character_sum = sum(1 - 2 * parity[u & x] for x in support)
            if half != zero_term - character_sum:
                raise AssertionError((n, truth, u, "different direct Walsh value"))
        sign = global_sign * (1 - 2 * parity[u & xor_support])
        normalized = sign * half
        normalized_corrected = sign * (half - zero_term)
        if not (normalized % 4 == normalized_corrected % 4 == 1):
            raise AssertionError((n, truth, u, "affine sign normalization failed"))
        raw_q = ((normalized - 1) // 4) % 2
        magnitude_q = int(abs(half) % 8 in (3, 5))
        if raw_q != magnitude_q:
            raise AssertionError((n, truth, u, "magnitude residue failed"))
        corrected_q = ((normalized_corrected - 1) // 4) % 2
        expected_q = constant ^ parity[u & linear_mask]
        for pair in quadratic_pairs:
            expected_q ^= int(u & pair == pair)
        if corrected_q != expected_q or (raw_q ^ corrected_q) != (n == 3 and u == 0):
            raise AssertionError((n, truth, u, "quadratic formula failed"))
        raw_values.append(raw_q)
        corrected_values.append(corrected_q)
        quadratic_pattern |= corrected_q << u
        # Digest actual signed and unnormalized integer Walsh/2 coefficients.
        digest.update(half.to_bytes(2, "little", signed=True))
    for name, values in (("raw", raw_values), ("corrected", corrected_values)):
        degree = algebraic_degree(values, axes, size)
        expected_max = 3 if name == "raw" and n == 3 else 2
        if degree > expected_max or (expected_max == 3 and degree != 3):
            raise AssertionError((n, truth, name, degree, "degree mismatch"))
        audit[name][degree] += 1
    return quadratic_pattern


def run_dimension(n, exhaustive, count, rng):
    import hashlib

    size = 1 << n
    parity = parity_table(size)
    axes = [1 << i for i in range(n)]
    pairs = [i | j for i, j in combinations(axes, 2)]
    audit = {"raw": Counter(), "corrected": Counter()}
    patterns = set()
    digest = hashlib.sha256()
    if exhaustive:
        # Set the last truth-table bit exactly when the other bits have even weight.
        truths = (tail | ((1 ^ (bin(tail).count("1") % 2)) << (size - 1)) for tail in range(1 << (size - 1)))
        count = 1 << (size - 1)
    else:
        tails = (rng.getrandbits(size - 1) for _ in range(count))
        truths = (tail | ((1 ^ (bin(tail).count("1") % 2)) << (size - 1)) for tail in tails)
    checked_direct = count if exhaustive else min(count, 8)
    for index, truth in enumerate(truths):
        pattern = check_table(n, truth, parity, axes, pairs, index < checked_direct, audit, digest)
        patterns.add(pattern)
    record = {
        "dimension": n,
        "mode": "exhaustive" if exhaustive else "deterministic seeded sample",
        "tables": count,
        "frequencies_checked": count * size,
        "direct_sum_tables": checked_direct,
        "direct_sum_frequencies": checked_direct * size,
        "distinct_corrected_quadratics": len(patterns),
        "raw_degree_histogram": dict(audit["raw"]),
        "corrected_degree_histogram": dict(audit["corrected"]),
        "true_walsh_half_sha256": digest.hexdigest(),
    }
    if exhaustive:
        distances = Counter()
        for left, right in combinations(sorted(patterns), 2):
            distances[bin(left ^ right).count("1")] += 1
        if any(distance < size // 4 for distance in distances):
            raise AssertionError((n, "distance smaller than RM(2,n) minimum"))
        record["distinct_quadratic_pair_distance_histogram"] = dict(distances)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("exhaustive", "sample"))
    args = parser.parse_args()
    rng = random.Random(SEED)
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    if args.mode == "exhaustive":
        records = [run_dimension(n, True, 0, rng) for n in EXHAUSTIVE_DIMENSIONS]
    else:
        records = [run_dimension(n, False, count, rng) for n, count in SAMPLE_COUNTS.items()]
    cpu_seconds = time.process_time() - cpu_start
    wall_seconds = time.perf_counter() - wall_start
    print(json.dumps({
        "status": "pass",
        "mode": args.mode,
        "seed": None if args.mode == "exhaustive" else SEED,
        "arithmetic": "exact Python integers",
        "python": platform.python_version(),
        "host": platform.node(),
        "machine": platform.machine(),
        "platform": platform.platform(),
        "wall_seconds": wall_seconds,
        "process_cpu_seconds": cpu_seconds,
        "mean_cpu_percent_of_one_core": 100 * cpu_seconds / wall_seconds,
        "records": records,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
