"""Integer verification of four small random affine basis examples."""

import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, product


def binary_rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            leading_bit = vector.bit_length()
            if leading_bit in pivots:
                vector ^= pivots[leading_bit]
            else:
                pivots[leading_bit] = vector
                break
    return len(pivots)


def all_nonempty_subsets(size):
    for cardinality in range(1, size + 1):
        yield from combinations(range(size), cardinality)


def basis_audit(bands, allocation_filter=lambda counts: True):
    entries = [(i, point) for i, band in enumerate(bands) for point in band]
    rank = binary_rank(point for _, point in entries)
    allocations = Counter()
    selected = 0
    appearances = [0] * len(entries)
    for basis in combinations(range(len(entries)), rank):
        if binary_rank(entries[index][1] for index in basis) != rank:
            continue
        counts = [0] * len(bands)
        for index in basis:
            counts[entries[index][0]] += 1
        allocations[tuple(counts)] += 1
        if allocation_filter(counts):
            selected += 1
            for index in basis:
                appearances[index] += 1
    marginals = [Fraction(appearance, selected) for appearance in appearances]
    summary = {
        "size": len(entries),
        "rank": rank,
        "bases": sum(allocations.values()),
        "allocations": len(allocations),
        "selected_bases": selected,
        "marginals": [str(marginal) for marginal in marginals],
        "squared_cost": str(sum(marginal**2 for marginal in marginals)),
    }
    return summary, set(allocations)


def predicted_allocations(bands, rank, bound):
    bounds = [
        (indices, bound(indices))
        for indices in all_nonempty_subsets(len(bands))
    ]
    return {
        counts
        for counts in product(*(range(len(band) + 1) for band in bands))
        if sum(counts) == rank
        and all(sum(counts[i] for i in indices) <= limit for indices, limit in bounds)
    }


def walsh_support(values):
    return [
        frequency
        for frequency in range(len(values))
        if sum(
            value * (-1) ** ((index & frequency).bit_count() % 2)
            for index, value in enumerate(values)
        )
    ]


def main():
    dependent = [[8 + 16 * i + u for u in range(4)] for i in range(4)]
    independent_directions = [list(range(4))] * 2 + [[0, 1, 4, 5]] * 2
    independent = [
        [2 ** (3 + i) + u for u in directions]
        for i, directions in enumerate(independent_directions)
    ]
    dependent_summary, dependent_allocations = basis_audit(dependent)
    independent_summary, independent_allocations = basis_audit(independent)
    for summary, expected_rank in ((dependent_summary, 5), (independent_summary, 7)):
        assert summary["rank"] == expected_rank
        assert summary["marginals"] == [str(Fraction(expected_rank, 16))] * 16

    assert dependent_allocations == predicted_allocations(
        dependent,
        5,
        lambda indices: 2 + binary_rank(dependent[i][0] for i in indices),
    )
    assert independent_allocations == predicted_allocations(
        independent,
        7,
        lambda indices: len(indices)
        + binary_rank(u for i in indices for u in independent_directions[i]),
    )

    offset_examples = {
        "parallel_nonuniform": [*range(8, 15), 24],
        "boolean_nonuniform": [*range(8, 16), 24],
    }
    results = {
        "python": sys.version.split()[0],
        "independent": independent_summary,
        "dependent": dependent_summary,
    }
    for name, offsets in offset_examples.items():
        bands = [[offset, offset ^ 32] for offset in offsets]
        summary, allocations = basis_audit(bands, lambda counts: counts[-1] == 1)
        assert summary["rank"] == 6
        assert allocations == predicted_allocations(
            bands, 6, lambda indices: 1 + binary_rank(offsets[i] for i in indices)
        )
        expected = [str(Fraction(5, 2 * (len(bands) - 1)))] * (2 * (len(bands) - 1))
        assert summary["marginals"] == expected + ["1/2"] * 2
        results[name] = summary

    truth_table = [
        (-1 if data % 8 == 0 else 1) if address == 0 else (-1) ** (data >> 3)
        for address in range(2)
        for data in range(16)
    ]
    support = walsh_support(truth_table)
    assert support == [*range(9), *range(16, 25)]
    results["boolean_nonuniform"]["walsh_support"] = support
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
