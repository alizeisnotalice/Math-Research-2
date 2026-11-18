#!/usr/bin/env python3
"""Exact Walsh examples and fractional-basis costs for quotient separators."""

from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
import time


def basis_add(basis, value):
    value = reduce_binary(basis, value)
    if value:
        basis[value.bit_length()] = value


def reduce_binary(basis, value):
    for bit_length in sorted(basis, reverse=True):
        if value & (1 << (bit_length - 1)):
            value ^= basis[bit_length]
    return value


def rank_binary(vectors):
    basis = {}
    for value in vectors:
        basis_add(basis, value)
    return len(basis)


def linear_subspaces(n):
    subspaces = {frozenset({0})}
    for vector in range(1, 1 << n):
        for subspace in tuple(subspaces):
            if vector not in subspace:
                subspaces.add(subspace | {vector ^ point for point in subspace})
    return subspaces


def nonempty_subsets(values):
    for cardinality in range(1, len(values) + 1):
        yield from itertools.combinations(values, cardinality)


def fractional_basis_cost(values):
    """Construct an exact minimum-norm base marginal and check its certificate."""
    nonloops = tuple(index for index, value in enumerate(values) if value)
    subset_ranks = {
        frozenset(subset): rank_binary(values[index] for index in subset)
        for subset in nonempty_subsets(nonloops)
    }
    old = frozenset()
    old_rank = 0
    remaining = frozenset(nonloops)
    marginals = [Fraction(0)] * len(values)
    for_rank = rank_binary(values)
    while remaining:
        cost, _, chosen = min(
            (Fraction(subset_ranks[old | chosen] - old_rank, len(chosen)), -len(chosen), tuple(sorted(chosen)))
            for chosen in subset_ranks
            if chosen <= remaining
        )
        for index in chosen:
            marginals[index] = cost
        old |= frozenset(chosen)
        old_rank = subset_ranks[old]
        remaining -= frozenset(chosen)
    assert sum(marginals) == for_rank
    assert all(sum(marginals[index] for index in subset) <= rank for subset, rank in subset_ranks.items())
    norm = sum((value * value for value in marginals), start=Fraction(0))
    assert all(
        sum(marginals[index] for index in subset) >= norm
        for subset, rank in subset_ranks.items()
        if len(subset) == rank == for_rank
    )
    return norm, marginals


def quotient_report(support, dimension):
    support = tuple(support)
    summaries = {}
    for core in linear_subspaces(dimension):
        basis = {}
        for value in sorted(core):
            basis_add(basis, value)
        quotient_values = [reduce_binary(basis, value) for value in support]
        rank = rank_binary(quotient_values)
        touches = 0 in quotient_values
        cost, marginals = fractional_basis_cost(quotient_values)
        key = (touches, rank)
        row = summaries.setdefault(key, {"count": 0, "costs": {}, "maximum_example": None})
        row["count"] += 1
        row["costs"][str(cost)] = row["costs"].get(str(cost), 0) + 1
        current = row["maximum_example"]
        if current is None or cost > Fraction(current["cost"]):
            row["maximum_example"] = {
                "cost": str(cost), "core": sorted(core),
                "marginals": [str(value) for value in marginals],
            }
    return [dict(row, touches=touches, rank=rank) for (touches, rank), row in sorted(summaries.items())]


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


def boundary_support():
    values = []
    for point in range(16):
        w, v, y, z = ((point >> coordinate) & 1 for coordinate in range(4))
        values.append(1 - 2 * (w & (v ^ (y & z))))
    support = walsh(values)
    assert set(support) == {0, 1, 2, 3, 6, 7, 10, 11, 14, 15}
    return tuple(value | 16 for value in support)


def addressed_example(depths, touching):
    address_bits = max(depths)
    address_count = 1 << address_bits
    selectors = []
    next_address = 0
    for depth in depths:
        cell_size = 1 << (address_bits - depth)
        assert next_address % cell_size == 0
        selectors.extend([len(selectors)] * cell_size)
        next_address += cell_size
    assert len(selectors) == address_count
    cell_indicators = [tuple(selector == cell for selector in selectors) for cell in sorted(set(selectors))]
    assert len(cell_indicators) == len(touching)
    private_dimensions = [2 if touch else 3 for touch in touching]
    private_starts = [sum(private_dimensions[:index]) for index in range(len(touching))]
    data_bits = sum(private_dimensions)
    cells = {cell: index for index, cell in enumerate(sorted(set(selectors)))}

    def value(address, data):
        cell = cells[selectors[address]]
        bits = data >> private_starts[cell]
        exponent = (bits & 1) & ((bits >> 1) & 1)
        if not touching[cell]:
            exponent ^= (bits >> 2) & 1
        return 1 - 2 * exponent

    true_spectrum = walsh(
        value(address, data) for data in range(1 << data_bits) for address in range(address_count)
    )
    expected_spectrum = {}
    for index, indicator in enumerate(cell_indicators):
        for u_frequency, u_entry in walsh(indicator).items():
            for local_label in range(4):
                sign = -1 if local_label == 3 else 1
                local_label |= 0 if touching[index] else 4
                total_label = local_label << private_starts[index]
                coefficient = u_entry * (1 << (data_bits - 1)) * sign
                total_frequency = (total_label << address_bits) | u_frequency
                expected_spectrum[total_frequency] = expected_spectrum.get(total_frequency, 0) + coefficient
    assert true_spectrum == {frequency: value for frequency, value in expected_spectrum.items() if value}
    nonglobal_zero = {frequency >> address_bits for frequency in true_spectrum if frequency >> address_bits}
    groups = [
        {label for label in nonglobal_zero if label >> start and not label >> (start + bits)}
        for start, bits in zip(private_starts, private_dimensions)
    ]
    assert all(len(group) == (3 if touch else 4) for group, touch in zip(groups, touching))
    projected = {0} | nonglobal_zero
    edges = [four for four in itertools.combinations(sorted(projected), 4) if four[0] ^ four[1] ^ four[2] ^ four[3] == 0]
    assert all(any(set(edge) - {0} <= group for group in groups) for edge in edges)
    assert all(any(group <= set(edge) for edge in edges) for group in groups)
    group_rows = []
    seen = set()
    for index, group in enumerate(groups):
        probability = Fraction(sum(cell_indicators[index]), address_count)
        increment = rank_binary(seen | group) - rank_binary(seen)
        group_band_points = sum(frequency >> address_bits in group for frequency in true_spectrum)
        assert group_band_points == len(group) / probability
        cost = Fraction(increment**2, group_band_points) / probability
        assert cost == (Fraction(4, 3) if touching[index] else Fraction(9, 4))
        group_rows.append({"p": str(probability), "r": increment, "q": group_band_points, "squared_constant": str(cost)})
        seen.update(group)
    return {
        "depths": depths,
        "touching": touching,
        "truth_table_size": 1 << (data_bits + address_bits),
        "global_sparsity": len(true_spectrum),
        "zero_label_band_size": sum(frequency >> address_bits == 0 for frequency in true_spectrum),
        "separated_groups": group_rows,
    }


def main():
    started = time.monotonic()
    result = {
        "schema": "a28-r067-separator-budget-v1",
        "examples": [
            addressed_example((1, 2, 3, 3), (True, True, True, True)),
            addressed_example((1, 2, 3, 3), (False, False, False, False)),
            addressed_example((1, 2, 3, 3), (True, False, False, True)),
        ],
        "boundary_quotient_profiles": quotient_report(boundary_support(), 5),
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers and rational arithmetic; no randomness",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.monotonic() - started,
    }
    Path(sys.argv[1]).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
