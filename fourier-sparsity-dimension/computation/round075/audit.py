"""Exact simultaneous physical and collision-free compression for two Boolean branches."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations_with_replacement
import json
from pathlib import Path
import platform
import time


def bits(mask):
    while mask:
        first = mask & -mask
        yield first.bit_length() - 1
        mask ^= first


def walsh(signs):
    signs = list(signs)
    step = 1
    while step < len(signs):
        for base in range(0, len(signs), 2 * step):
            for i in range(base, base + step):
                left, right = signs[i], signs[i + step]
                signs[i], signs[i + step] = left + right, left - right
        step *= 2
    return signs


def span(vectors):
    space = {0}
    for v in vectors:
        space |= {x ^ v for x in space}
    return frozenset(space)


@lru_cache(maxsize=None)
def vector_rank(label_mask):
    return len(span(bits(label_mask))).bit_length() - 1


@lru_cache(maxsize=None)
def linear_spaces(n):
    new = {frozenset({0})}
    all_spaces = set(new)
    for _ in range(n):
        new = {space | {v ^ u for u in space} for space in new for v in range(1 << n) if v not in space}
        all_spaces |= new
    return sorted(all_spaces, key=lambda x: (len(x), sorted(x)))


@lru_cache(maxsize=None)
def affine_spaces(n):
    spaces = {tuple(sorted(v ^ u for u in space))
              for space in linear_spaces(n) for v in range(1 << n)}
    return sorted(spaces, key=lambda x: (len(x), x))


def label_support(values):
    return sum(1 << label for label, value in enumerate(walsh(values)) if value)


def combine_masks(masks):
    union = 0
    for mask in masks:
        union |= mask
    return union


def collision_budget(transform, query_size):
    basis = frozenset({0})
    differences = 0
    multiplicity = 0
    for label in range(1, len(transform) // query_size):
        row = [j for j in range(query_size) if transform[j + label * query_size]]
        if not row or label in basis:
            continue
        basis |= {v ^ label for v in basis}
        multiplicity += len(row) - 1
        for frequency in row[1:]:
            differences |= 1 << (frequency ^ row[0])
    return differences, multiplicity


def truth_rank_profiles(n, k):
    size = 1 << n
    query_size = 1 << k
    data_size = 1 << (n - k)
    P = affine_spaces(k)
    dims = [len(w).bit_length() - 1 for w in P]
    by_dim = [sum(1 << i for i, dim in enumerate(dims) if dim == h) for h in range(k + 1)]
    actual = Counter()
    actual_example = {}
    avoid = Counter()
    avoid_example = {}

    @lru_cache(maxsize=None)
    def local(mask):
        return label_support(1 - 2 * ((mask >> i) & 1) for i in range(data_size))

    for mask in range(1 << size):
        truth = [1 - 2 * ((mask >> i) & 1) for i in range(size)]
        transform = walsh(truth)
        address_masks = []
        for row in range(query_size):
            row_mask = sum(((mask >> (i * query_size + row)) & 1) << i for i in range(data_size))
            address_masks.append(local(row_mask))
        whole_support = combine_masks(address_masks)
        assert whole_support == sum(1 << label for label in range(data_size)
                                   if any(transform[label * query_size + i] for i in range(query_size)))
        r = vector_rank(whole_support)
        good = sum(1 << i for i, w in enumerate(P)
                   if vector_rank(combine_masks(address_masks[j] for j in w)) == r)
        assert good & (1 << (len(P) - 1))
        d, budget = collision_budget(transform, query_size)
        actual[(good, budget)] += 1
        actual_example.setdefault((good, budget), mask)
        avoid[d] += 1
        avoid_example.setdefault(d, mask)

    @lru_cache(maxsize=None)
    def physical_min(good):
        return dims[(good & -good).bit_length() - 1]

    @lru_cache(maxsize=None)
    def avoidance_min(diffs):
        for W in reversed(linear_spaces(k)):
            if not any((diffs >> w) & 1 for w in W):
                return k + 1 - len(W).bit_length()
        raise AssertionError("zero kernel always qualifies")

    actual_distribution = Counter()
    actual_max = {}
    nonzero_probability_equalities = 0
    for first, second in combinations_with_replacement(sorted(actual), 2):
        good0, budget0 = first
        good1, budget1 = second
        common = good0 & good1
        h0, h1 = sorted([physical_min(good0), physical_min(good1)])
        h = physical_min(common)
        assert max(h0, h1) <= h <= min(k, h0 + h1 + 1)
        pair_weight = actual[first] * actual[second] * (1 if first == second else 2)
        actual_distribution[(h0, h1, h)] += pair_weight
        for d, eligible in enumerate(by_dim):
            total = eligible.bit_count()
            failed = total - (common & eligible).bit_count()
            left = failed * 2 * (query_size - 1)
            right = total * (budget0 + budget1) * ((query_size >> d) - 1)
            assert left <= right
            if failed and left == right:
                nonzero_probability_equalities += pair_weight
        gap = h - h1
        if gap > actual_max.get("gap", -1):
            actual_max = {"independent_h": [h0, h1], "joint_h": h, "gap": gap,
                          "truth_masks": [actual_example[first], actual_example[second]]}

    avoidance_distribution = Counter()
    avoidance_max = {}
    for first, second in combinations_with_replacement(sorted(avoid), 2):
        h0, h1 = sorted([avoidance_min(first), avoidance_min(second)])
        union = first | second
        h = avoidance_min(union)
        assert max(h0, h1) <= h <= min(k, h0 + h1)
        assert h <= (union.bit_count() + 1).bit_length() - 1
        avoidance_distribution[(h0, h1, h)] += avoid[first] * avoid[second] * (1 if first == second else 2)
        gap = h - h1
        if gap > avoidance_max.get("gap", -1):
            avoidance_max = {"independent_h": [h0, h1], "joint_h": h, "gap": gap,
                             "truth_masks": [avoid_example[first], avoid_example[second]],
                             "difference_masks": [first, second]}

    expected_pairs = 1 << (2 * size)
    assert sum(actual_distribution.values()) == sum(avoidance_distribution.values()) == expected_pairs
    return {"dimension": n, "query_dimension": k, "all_ordered_branch_pairs": expected_pairs,
            "truth_tables": 1 << size, "affine_restrictions": len(P),
            "rank_retention_and_budget_profiles": len(actual), "difference_profiles": len(avoid),
            "positive_failure_probability_bound_equalities": nonzero_probability_equalities,
            "actual_simultaneous_worst_gap": actual_max,
            "collision_free_worst_gap": avoidance_max,
            "actual_h0_h1_joint_histogram": {",".join(map(str, key)): count
                                              for key, count in sorted(actual_distribution.items())},
            "collision_free_h0_h1_joint_histogram": {",".join(map(str, key)): count
                                                      for key, count in sorted(avoidance_distribution.items())}}


def bent(y, a):
    bit = sum(((y >> j) & 1) * ((y >> (j + 1)) & 1) for j in range(0, a, 2)) & 1
    return 1 - 2 * bit


def balanced_address_families(k):
    query_bits = 2 * k + 2
    size = 1 << query_bits
    m = 1 << k
    guard = 1 << (2 * k)
    A = [dict(zip(range(m), range(m))), {guard ^ (j << k): j for j in range(m)}]
    branch_data = 1 << m

    def leaf(b, y, z):
        if y in A[b]:
            return 1 - 2 * ((z >> A[b][y]) & 1)
        return bent(y, query_bits)

    core = [walsh(leaf(b, y, z) for z in range(branch_data) for y in range(size))
            for b in range(2)]
    for spectrum in core:
        assert all(spectrum[y] for y in range(size))
        assert sum(bool(c) for c in spectrum) == (1 + m) * size

    ext_size = 2 * branch_data ** 2

    def F(y, outer):
        branch = outer & 1
        return leaf(branch, y, (outer >> (1 + branch * m)) & (branch_data - 1))

    truth = [F(y, outer) for outer in range(ext_size) for y in range(size)]
    spectrum = walsh(truth)
    assert all(spectrum[y] for y in range(size))
    branch_rows = [[label_support(leaf(b, y, z) for z in range(branch_data)) for y in range(size)]
                   for b in range(2)]
    joint_rows = [label_support(F(y, x) for x in range(ext_size)) for y in range(size)]
    assert vector_rank(combine_masks(joint_rows)) == 1 + 2 * m
    assert vector_rank(combine_masks(branch_rows[0])) == vector_rank(combine_masks(branch_rows[1])) == m
    all_P = affine_spaces(query_bits)
    seq = [[], [], []]
    hist = Counter()
    for P in all_P:
        ranks = [vector_rank(combine_masks(rows[y] for y in P)) for rows in (*branch_rows, joint_rows)]
        counts = [sum(y in A[b] for y in P) for b in range(2)]
        assert ranks[:2] == counts
        assert ranks[2] == ((1 + sum(counts)) if any(counts) else 0)
        hist[(len(P).bit_length() - 1, *ranks)] += 1
        for target, rank in enumerate(ranks):
            if rank == (m if target < 2 else 1 + 2 * m):
                seq[target].append(P)
    min_dims = [len(rows[0]).bit_length() - 1 for rows in seq]
    assert min_dims == [k, k, 2 * k + 1]
    minimal_coaffine = seq[2][0]
    assert minimal_coaffine == tuple(range(1 << (2 * k + 1)))
    retained = walsh(F(y, x) for x in range(ext_size) for y in minimal_coaffine)
    actual_bands = {i // len(minimal_coaffine) for i, c in enumerate(retained) if c}
    assert vector_rank(sum(1 << i for i in actual_bands)) == 1 + 2 * m
    return {"k": k, "query_bits": query_bits, "all_query_flats": len(all_P),
            "active_query_flats": [sorted(row) for row in A],
            "branch_tail_counts": [m * size] * 2, "global_full_truth_table_entries": len(truth),
            "joint_zero_band_size": size,
            "joint_tail_count": sum(bool(c) for c in spectrum[size:]),
            "full_exterior_rank": 1 + 2 * m,
            "branch_and_joint_minimum_compressed_query_dims": min_dims,
            "joint_witness_query_flat": list(minimal_coaffine),
            "joint_witness_exact_sparse_restriction": sum(bool(c) for c in retained),
            "combined_rank_histogram_h_r0_r1_r": {",".join(map(str,key)): count
                                                  for key, count in sorted(hist.items())}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    start = time.perf_counter()
    cases = [(n, k) for n in range(1, 5) for k in range(n + 1)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        pending = [pool.submit(truth_rank_profiles, n, k) for n, k in cases]
        records = [process.result() for process in pending]
    exact = [balanced_address_families(k) for k in range(3)]
    result = {"schema": "a28-r074-simultaneous-compression-v1", "python": platform.python_version(),
              "workers": args.workers, "finite_truth_census": records,
              "sharp_balanced_construction": exact, "elapsed_seconds": time.perf_counter() - start}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    for x in records:
        print(x["dimension"], x["query_dimension"], x["actual_simultaneous_worst_gap"],
              x["collision_free_worst_gap"])
    print(result["elapsed_seconds"])
