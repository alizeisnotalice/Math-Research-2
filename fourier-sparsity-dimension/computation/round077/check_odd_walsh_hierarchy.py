#!/usr/bin/env python3
"""Exact pair census for the first differing residue of normalized Walsh/2."""

from collections import Counter
from itertools import combinations
import json
import platform
import time

from check_odd_walsh_mod8 import algebraic_degree, walsh_half


def verify_dimension(n):
    size = 1 << n
    axes = [1 << i for i in range(n)]
    odd_count = 1 << (size - 1)
    representatives = []
    for truth in range(1 << size):
        if bin(truth).count("1") % 4 != 3:
            continue
        support_xor = 0
        for x in range(size):
            if truth >> x & 1:
                support_xor ^= x
        if support_xor != 0:
            continue
        spectrum = walsh_half(size, truth)
        if any(value % 4 != 1 for value in spectrum):
            raise AssertionError((n, truth, "unexpected mod-4 Walsh normalization"))
        planes = tuple(sum(((value >> bit) & 1) << u for u, value in enumerate(spectrum))
                       for bit in range(2, n))
        representatives.append(planes)
    if len(representatives) != odd_count // (2 * size):
        raise AssertionError((n, len(representatives), "incomplete normalized orbit representatives"))
    if len(set(representatives)) != len(representatives):
        raise AssertionError((n, "duplicate normalized spectra"))
    first_split_counts = Counter()
    degree_cache = {}
    for left, right in combinations(representatives, 2):
        for bit, (left_plane, right_plane) in enumerate(zip(left, right), start=2):
            pattern = left_plane ^ right_plane
            if pattern == 0:
                continue
            if pattern not in degree_cache:
                values = [(pattern >> u) & 1 for u in range(size)]
                degree_cache[pattern] = algebraic_degree(values, axes, size)
            degree = degree_cache[pattern]
            distance = bin(pattern).count("1")
            if degree > bit or distance < size >> bit:
                raise AssertionError((n, bit, pattern, degree, distance, "unexpected first split"))
            first_split_counts[bit, degree, distance] += 1
            break
        else:
            raise AssertionError((n, left, right, "distinct representatives never split"))
    return {
        "dimension": n,
        "normalization_orbit_representatives": len(representatives),
        "distinct_pairs": len(representatives) * (len(representatives) - 1) // 2,
        "first_split": [dict(bit=bit, degree=degree, distance=distance, pairs=pairs)
                        for (bit, degree, distance), pairs in sorted(first_split_counts.items())],
    }


def main():
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    records = [verify_dimension(n) for n in (3, 4)]
    print(json.dumps({
        "status": "pass",
        "python": platform.python_version(),
        "host": platform.node(),
        "arithmetic": "exact Python integers",
        "records": records,
        "wall_seconds": time.perf_counter() - wall_start,
        "process_cpu_seconds": time.process_time() - cpu_start,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
