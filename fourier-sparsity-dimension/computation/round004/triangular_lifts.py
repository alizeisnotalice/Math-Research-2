#!/usr/bin/env python3
"""Exact sharp sign-code lifts for address spectra; no numerical tolerance."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import platform
import random
import sys


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p in pivots:
                row ^= pivots[p]
            else:
                pivots[p] = row
                break
    return len(pivots)


def walsh(v):
    v = v[:]
    step = 1
    while step < len(v):
        for start in range(0, len(v), 2 * step):
            for j in range(start, start + step):
                a, b = v[j], v[j + step]
                v[j], v[j + step] = a + b, a - b
        step *= 2
    return v


def check(k):
    t = 1 << k
    # Variables are t intercepts eta_i and sum_j 2^j triangular parameters.
    slots = [(j, prefix) for j in range(k) for prefix in range(1 << j)]
    assert len(slots) == t - 1
    support, coeffs = [], []
    for i in range(t):
        for a in range(t):
            frequency = 1 << i
            for slot, (j, prefix) in enumerate(slots):
                if (i & ((1 << j) - 1)) == prefix and ((a >> j) & 1):
                    frequency |= 1 << (t + slot)
            support.append(frequency)
            coeffs.append(-1 if (a & i).bit_count() & 1 else 1)
    d = 2 * t - 1
    assert len(set(support)) == t * t
    assert rank(support) == d
    assert rank([v | (1 << d) for v in support]) == d
    # All nonzero lifted differences cancel, an exact complete Boolean certificate.
    pair_sums = defaultdict(int)
    for i, a in enumerate(support):
        for j in range(i):
            pair_sums[a ^ support[j]] += coeffs[i] * coeffs[j]
    assert all(value == 0 for value in pair_sums.values())
    assert sum(v * v for v in coeffs) == t * t
    direct_truth = d <= 15
    if direct_truth:
        spectrum = [0] * (1 << d)
        for a, c in zip(support, coeffs):
            spectrum[a] = c
        assert all(abs(v) == t for v in walsh(spectrum))
    parameter_count = t - 1
    rng = random.Random(20261116 + k)
    params = range(1 << parameter_count) if k <= 3 else [rng.getrandbits(parameter_count) for _ in range(2000)]
    count = 0
    for z in params:
        image = []
        for i in range(t):
            b = i
            for slot, (j, prefix) in enumerate(slots):
                if (i & ((1 << j) - 1)) == prefix and ((z >> slot) & 1):
                    b ^= 1 << j
            image.append(b)
        assert len(set(image)) == t
        count += 1
    return {"k": k, "s": t * t, "linear_rank": d, "affine_rank": d - 1,
            "signcode_dimension": d, "pair_difference_classes": len(pair_sums),
            "all_pair_cancellations_checked": True, "full_inverse_walsh": direct_truth,
            "permutations_checked": count, "permutation_exhaustive": k <= 3,
            "spectrum_sha256": hashlib.sha256(json.dumps(list(zip(support, coeffs))).encode()).hexdigest()}


parser = argparse.ArgumentParser()
parser.add_argument("--max-k", type=int, default=5)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
result = {"method": "exact_F2_rank_all_pair_cancellations_and_triangular_permutations",
          "python": sys.version, "platform": platform.platform(), "tolerance": 0,
          "seed": 20261116, "rows": [check(k) for k in range(args.max_k + 1)],
          "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
args.out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "passed", "max_k": args.max_k, "largest_support": result['rows'][-1]['s']}))
