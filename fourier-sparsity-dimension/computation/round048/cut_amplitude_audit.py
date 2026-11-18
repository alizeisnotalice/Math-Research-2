#!/usr/bin/env python3
"""Exact parent-Fourier cut census for sparse selector families."""

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


C_VALUES = (Fraction(3, 2), Fraction(2), Fraction(5, 2), Fraction(3))


def parity(value):
    return int(value).bit_count() & 1


def walsh(values):
    result = list(values)
    width = 1
    while width < len(result):
        for start in range(0, len(result), 2 * width):
            for index in range(start, start + width):
                left, right = result[index], result[index + width]
                result[index], result[index + width] = left + right, left - right
        width *= 2
    return result


def rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


def path_band(node):
    numbers = {0: 1}
    while node:
        parent = (node - 1) // 2
        branch_sign = 1 if node == 2 * parent + 1 else -1
        numbers.update({frequency | (1 << parent): coefficient * branch_sign
                        for frequency, coefficient in numbers.items()})
        node = parent
    return numbers


def independent_tree(height, replacement_depth=None):
    bottom = (1 << (height - 1)) - 1
    original_nodes = (1 << height) - 1
    scale = 1 << (height - 1)
    target = None
    if replacement_depth is not None:
        target = (1 << replacement_depth) - 1
        scale = 1 << max(height - 1, replacement_depth + 4)

    def above(ancestor, node):
        while node > ancestor:
            node = (node - 1) // 2
        return node == ancestor

    spectrum = {
        frequency | (1 << node): sign * (scale >> (height - 1))
        for node in range(bottom, original_nodes)
        if target is None or not above(target, node)
        for frequency, sign in path_band(node).items()
    }
    mask = 0x0117
    if target is not None:
        local = walsh(1 - 2 * ((mask >> point) & 1) for point in range(16))
        for address, sign in path_band(target).items():
            for beta, coefficient in enumerate(local):
                if coefficient:
                    assert address | (beta << original_nodes) not in spectrum
                    spectrum[address | (beta << original_nodes)] = sign * coefficient * (scale >> (4 + replacement_depth))

    def value(point):
        node = 0
        while node < bottom:
            if node == target:
                return 1 - 2 * ((mask >> (point >> original_nodes)) & 1)
            node = 2 * node + 1 + ((point >> node) & 1)
        return 1 - 2 * ((point >> node) & 1)

    return original_nodes + (4 if target is not None else 0), spectrum, scale, value


def addressing(address_bits):
    count = 1 << address_bits
    spectrum = {
        frequency | (1 << (address_bits + label)): 1 - 2 * parity(frequency & label)
        for label in range(count)
        for frequency in range(count)
    }

    def value(point):
        return 1 - 2 * ((point >> (address_bits + (point & (count - 1)))) & 1)

    return address_bits + count, spectrum, count, value


def child_spectra(spectrum, query):
    pivot = 1 << (query.bit_length() - 1)
    pairs = {}
    for frequency, coefficient in spectrum.items():
        side = int(bool(frequency & pivot))
        representative = frequency ^ (query if side else 0)
        pairs.setdefault(representative, [0, 0])[side] = coefficient
    left = {}
    right = {}
    overlap = 0
    equal_child = 0
    both_parent = 0
    for frequency, (first, second) in pairs.items():
        positive = first + second
        negative = first - second
        if positive:
            left[frequency] = positive
        if negative:
            right[frequency] = negative
        if positive and negative:
            overlap += 1
            equal_child += abs(positive) == abs(negative)
        both_parent += bool(first and second)
    assert len(spectrum) == 2 * (len(left) + len(right)) - 2 * overlap - equal_child
    return left, right, overlap, equal_child, both_parent


def payable(parent, left, right, intersection, constant):
    """Exactly compare sqrt(left)+sqrt(right) <= sqrt(parent)+intersection/constant."""
    shift = Fraction(intersection) / constant
    if left == right:
        remaining = 4 * left - parent - shift * shift
        return remaining <= 0 or remaining * remaining <= 4 * shift * shift * parent
    bias = parent + shift * shift - left - right
    if bias < 0 and 4 * shift * shift * parent <= bias * bias:
        return False
    remaining = 4 * left * right - bias * bias - 4 * shift * shift * parent
    square_limit = 16 * bias * bias * shift * shift * parent
    if bias >= 0:
        return remaining <= 0 or remaining * remaining <= square_limit
    return remaining < 0 and remaining * remaining >= square_limit


def direct_cut(dimension, scale, value, query, left, right):
    pivot = 1 << (query.bit_length() - 1)
    lower = pivot - 1
    representatives = [(point & lower) | ((point & ~lower) << 1)
                       for point in range(1 << (dimension - 1))]
    for bit, branch in enumerate((left, right)):
        transformed = walsh(value(point | pivot * (parity(point & query) ^ bit))
                            for point in representatives)
        assert all(coefficient * scale == branch.get(frequency, 0) * len(representatives)
                   for coefficient, frequency in zip(transformed, representatives))


def audit_model(name, model, check_c2):
    dimension, spectrum, scale, value = model
    size = len(spectrum)
    flag = 1 << dimension
    parent_augmented = rank(frequency | flag for frequency in spectrum)
    candidates = sorted({first ^ second for first, second in combinations(spectrum, 2)})
    flat = len(set(map(abs, spectrum.values()))) == 1
    if dimension <= 12:
        direct = walsh(value(point) for point in range(1 << dimension))
        assert all(coefficient * scale == spectrum.get(frequency, 0) * len(direct)
                   for frequency, coefficient in enumerate(direct))
        direct_queries = set(candidates[::max(1, len(candidates) // 6)]) | {1}
    else:
        direct_queries = set()

    hist = Counter()
    paid = {str(constant): [] for constant in C_VALUES}
    root = None
    for query in candidates:
        left, right, overlap, equal_child, double = child_spectra(spectrum, query)
        common_rank = (rank(frequency | flag for frequency in left) +
                       rank(frequency | flag for frequency in right) - parent_augmented + 1)
        assert common_rank >= 0
        if flat:
            assert overlap == equal_child == size - 2 * double
            assert len(left) == len(right) == size - 3 * double // 2
            assert double % 2 == 0
        if query in direct_queries:
            union_rank = rank(frequency | flag for frequency in left.keys() | right.keys())
            assert union_rank == parent_augmented - 1
            direct_cut(dimension, scale, value, query, left, right)
        sizes = sorted((len(left), len(right)))
        profile = (*sizes, overlap, equal_child, common_rank, double)
        hist[profile] += 1
        for constant in C_VALUES:
            if payable(size, *sizes, common_rank, constant):
                paid[str(constant)].append(query)
        if query == 1:
            root = dict(zip(("s0", "s1", "M", "equal_child", "h", "doubles"), profile))
    assert root is not None and 1 in paid["2"]
    if check_c2 is not None:
        assert set(paid["2"]) == check_c2(candidates)
    return {
        "name": name,
        "ambient_dimension": dimension,
        "sparsity": size,
        "affine_rank": parent_augmented - 1,
        "flat": flat,
        "strict_cut_count": len(candidates),
        "direct_cut_transform_checks": len(direct_queries),
        "counts_payable_by_C": {constant: len(queries) for constant, queries in paid.items()},
        "payable_C2_masks": paid["2"] if len(paid["2"]) <= 128 else paid["2"][:128],
        "root_cut": root,
        "metric_profile_order": ["min_child_s", "max_child_s", "M", "equal_child", "intersection_augmented_rank", "doubles"],
        "profile_counts": {",".join(map(str, profile)): count for profile, count in sorted(hist.items())},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    models = []
    for height in range(2, 6):
        expected = set if height == 2 else lambda _: {1}
        models.append(audit_model(f"independent_tree_h{height}", independent_tree(height), expected))
    for address_bits in range(1, 5):
        expected = set if address_bits <= 2 else lambda _, bits=address_bits: set(range(1, 1 << bits))
        models.append(audit_model(f"address_bits{address_bits}", addressing(address_bits), expected))
    for height, depth in ((3, 1), (4, 2), (5, 2)):
        models.append(audit_model(f"rare_leaf_h{height}_depth{depth}",
                                  independent_tree(height, depth), None))
    result = {
        "schema": "a28-r048-cut-amplitude-audit-v1",
        "models": models,
        "total_strict_cuts": sum(row["strict_cut_count"] for row in models),
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers and fractions; no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"checked {result['total_strict_cuts']} strict cuts in {len(models)} models")


if __name__ == "__main__":
    main()
