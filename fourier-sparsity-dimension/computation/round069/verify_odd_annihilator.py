#!/usr/bin/env python3
"""Exact checks for five-level means and the rank of pairing defects."""

from collections import Counter
from functools import reduce
import hashlib
import itertools
import json
from operator import xor
from pathlib import Path
import platform
import sys
import time


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for left in range(start, start + width):
                right = left + width
                first, second = values[left], values[right]
                values[left], values[right] = first + second, first - second
        width *= 2
    return {index: value for index, value in enumerate(values) if value}


def rank_binary(vectors):
    basis = {}
    for value in vectors:
        while value and value.bit_length() in basis:
            value ^= basis[value.bit_length()]
        if value:
            basis[value.bit_length()] = value
    return len(basis)


def linear_subspaces(n):
    subspaces = {frozenset({0})}
    for vector in range(1, 1 << n):
        for subspace in tuple(subspaces):
            if vector not in subspace:
                subspaces.add(subspace | {vector ^ point for point in subspace})
    return subspaces


def zero_sum_subsets(support, max_size):
    support = tuple(support)
    return [
        subset
        for cardinality in range(4, max_size + 1, 2)
        for subset in itertools.combinations(support, cardinality)
        if reduce(xor, subset, 0) == 0
    ]


def pair_matchings(labels):
    if not labels:
        yield ()
        return
    first = labels[0]
    for index, second in enumerate(labels[1:], 1):
        rest = labels[1:index] + labels[index + 1:]
        for matched in pair_matchings(rest):
            yield ((first, second), *matched)


def values_for_mean(numerator):
    if not numerator:
        return (1, -1, 1, -1)
    sign = 1 if numerator > 0 else -1
    if abs(numerator) == 4:
        return (sign,) * 4
    assert abs(numerator) == 2
    return (-sign, sign, sign, sign)


def ten_point_witness():
    d = (1, 2, 4, 8, 15)
    e = tuple(16 << i for i in range(5))
    phase = (1, 0, 0, 0, 0)
    labels = tuple(sorted(label for first, difference in zip(e, d) for label in (first, first ^ difference)))
    assert rank_binary(labels) == 9
    relations = zero_sum_subsets(labels, 10)
    assert relations == [labels]
    all_matchings = Counter(rank_binary(first ^ second for first, second in matching) for matching in pair_matchings(labels))
    assert all_matchings == {4: 945}

    activity = Counter()
    values = []
    for common in range(1 << 9):
        active = [index for index, difference in enumerate(d) if (common & difference).bit_count() % 2 == phase[index]]
        assert len(active) in (0, 2, 4)
        activity[len(active)] += 1
        values.append(sum(1 - 2 * ((common & e[index]) != 0) for index in active))
    mean_walsh = walsh(values)
    assert tuple(sorted(mean_walsh)) == labels
    assert all(abs(value) == 256 for value in mean_walsh.values())
    assert all(value * (value**2 - 4) * (value**2 - 16) == 0 for value in values)
    assert sum(abs(value) for value in mean_walsh.values()) == 2560
    truth_rows = [values_for_mean(numerator) for numerator in values]
    assert all(sum(row) == numerator for row, numerator in zip(truth_rows, values))
    private_slice_sparsities = Counter(len(walsh(row)) for row in truth_rows)
    assert set(private_slice_sparsities) == {1, 4}
    full_spectrum = walsh(row[private] for private in range(4) for row in truth_rows)
    assert {frequency: value for frequency, value in full_spectrum.items() if frequency < 1 << 9} == mean_walsh
    return {
        "common_phase": phase,
        "common_support": labels,
        "common_support_rank": 9,
        "local_relation_rank": 1,
        "perfect_matching_difference_rank_counts": all_matchings,
        "minimum_pairing_defect_rank": 4,
        "matching_chosen_by_pairs_on_e_i": list(d),
        "common_class_activity_by_query": activity,
        "four_times_common_mean_histogram": Counter(values),
        "private_slice_sparsity_histogram": private_slice_sparsities,
        "true_boolean_truth_table_size": len(truth_rows) * 4,
        "full_sparsity": len(full_spectrum),
        "full_rank": rank_binary(full_spectrum),
    }


def exceptional_absorption_witness():
    def value(address, exceptional, private):
        if address:
            return 1 - 2 * ((private & 1) & ((private >> 1) & 1))
        return 1 - 2 * (exceptional == 0)

    original_size = [
        len(walsh(value(address, exceptional, private) for private in range(4) for exceptional in range(8)))
        for address in range(2)
    ]
    promoted_size = Counter(
        len(walsh(value(address, exceptional, private) for private in range(4)))
        for address in range(2)
        for exceptional in range(8)
    )
    common_rows = [
        walsh(sum(value(address, exceptional, private) for private in range(4)) for exceptional in range(8))
        for address in range(2)
    ]
    assert set(common_rows[0]) == set(range(8)) and set(common_rows[1]) == {0}
    relations = zero_sum_subsets(common_rows[0], 10)
    acceptable_dimensions = Counter()
    for subspace in linear_subspaces(3):
        projected = [min(value ^ direction for direction in subspace) for value in range(8)]
        if all(all(count % 2 == 0 for count in Counter(projected[label] for label in r).values()) for r in relations):
            acceptable_dimensions[len(subspace).bit_length() - 1] += 1
    assert original_size == [8, 4]
    assert promoted_size == {1: 8, 4: 8}
    assert min(acceptable_dimensions) == 2
    return {
        "original_slice_sparsities": original_size,
        "promoted_exceptional_core_rank": 3,
        "promoted_slice_sparsity_histogram": promoted_size,
        "short_relations": len(relations),
        "permitted_defect_ranks": acceptable_dimensions,
        "minimum_pairing_defect_rank": min(acceptable_dimensions),
    }


def main():
    started = time.monotonic()
    result = {
        "schema": "a28-r069-odd-annihilator-v1",
        "ten_point_witness": ten_point_witness(),
        "exceptional_absorption": exceptional_absorption_witness(),
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers; no randomness",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.monotonic() - started,
    }
    Path(sys.argv[1]).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
