#!/usr/bin/env python3
"""Exact Walsh and local base costs for hierarchical activation examples."""

from fractions import Fraction
import hashlib
import itertools
import json
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


def rank(vectors):
    basis = {}
    for vector in vectors:
        while vector and vector.bit_length() in basis:
            vector ^= basis[vector.bit_length()]
        if vector:
            basis[vector.bit_length()] = vector
    return len(basis)


def affine_rank(points, dimension):
    return rank(point | (1 << dimension) for point in points) - 1


def full_spectrum(value, address_bits, data_bits):
    return walsh(
        value(address, data)
        for data in range(1 << data_bits)
        for address in range(1 << address_bits)
    )


def local_bands(value, address_bits, data_bits):
    results = {}
    for address in range(1 << address_bits):
        for label, amplitude in walsh(value(address, data) for data in range(1 << data_bits)).items():
            results.setdefault(label, {})[address] = amplitude
    return results


def band_checks(spectrum, bands, address_bits, laminar=True):
    address_count = 1 << address_bits
    reconstructed = {
        address_freq | (label << address_bits): value
        for label, activity in bands.items()
        for address_freq, value in walsh(activity.get(address, 0) for address in range(address_count)).items()
    }
    assert reconstructed == spectrum
    assert all(len(set(map(abs, activity.values()))) == 1 for activity in bands.values())
    assert all(
        len(walsh(activity.get(address, 0) for address in range(address_count))) * len(activity) == address_count
        for activity in bands.values()
    )
    domains = {frozenset(activity) for activity in bands.values()}
    assert all(not (left & right) or left <= right or right <= left for left, right in itertools.combinations(domains, 2)) == laminar
    terminal = {domain for domain in domains if not any(other < domain for other in domains)}
    assert set().union(*terminal) == set(range(address_count))
    assert sum(map(len, terminal)) == address_count
    assert len({len(block) for block in terminal}) == 1
    size = len(next(iter(terminal)))
    assert min(map(len, domains)) == size
    packing = Fraction(address_count, size)
    assert packing == len(terminal)
    return domains, terminal, packing


def independent_tree(height, bent):
    leaves = 1 << height
    address_bits = leaves - 1
    label_width = 3 if bent else 1
    data_bits = label_width * leaves

    def value(address, data):
        node = 1
        while node < leaves:
            node = 2 * node + ((address >> (node - 1)) & 1)
        block = data >> (label_width * (node - leaves))
        exponent = block & 1
        if bent:
            exponent ^= ((block >> 1) & 1) & ((block >> 2) & 1)
        return 1 - 2 * exponent

    spectrum = full_spectrum(value, address_bits, data_bits)
    bands = local_bands(value, address_bits, data_bits)
    domains, terminal, packing = band_checks(spectrum, bands, address_bits)
    assert len(domains) == len(terminal) == leaves
    assert len(spectrum) == (4 if bent else 1) * leaves**2
    assert packing == leaves
    augmented_plus_one = affine_rank(spectrum, address_bits + data_bits) + 2
    assert augmented_plus_one**2 == 4 * len(spectrum)
    return {
        "height": height,
        "leaf_function": "private_character_times_two_bit_bent" if bent else "private_character",
        "true_sparsity": len(spectrum),
        "linear_rank": rank(spectrum),
        "affine_rank_plus_two": augmented_plus_one,
        "activity_packing": str(packing),
    }


def nested_common_addresses(height):
    address_bits = height
    data_bits = (1 << (height + 1)) - 2

    def path_mask(address):
        labels = [((1 << depth) - 2) + (address & ((1 << depth) - 1)) for depth in range(1, height + 1)]
        return sum(1 << label for label in labels)

    paths = [path_mask(address) for address in range(1 << height)]

    def value(address, data):
        return 1 - 2 * (not (data & paths[address]))

    spectrum = full_spectrum(value, address_bits, data_bits)
    bands = local_bands(value, address_bits, data_bits)
    domains, terminal, packing = band_checks(spectrum, bands, address_bits)
    assert len(domains) == (1 << (height + 1)) - 1
    assert len(terminal) == 1 << height
    assert len(spectrum) == 1 + sum(1 << (3 * depth - 1) for depth in range(1, height + 1))
    assert rank(spectrum) == address_bits + data_bits
    assert all(
        sum(address in activity for activity in bands.values()) == 1 << height
        for address in range(1 << address_bits)
    )
    return {
        "height": height,
        "true_sparsity": len(spectrum),
        "linear_rank": rank(spectrum),
        "distinct_activity_domains": len(domains),
        "distinct_terminal_domains": len(terminal),
        "activity_packing": str(packing),
        "coactive_labels_per_query": 1 << height,
    }


def crossing_common_addresses(height):
    data_bits = 2 * height
    paths = [
        sum(1 << (2 * index + ((address >> index) & 1)) for index in range(height))
        for address in range(1 << height)
    ]

    def value(address, data):
        return 1 - 2 * (not (data & paths[address]))

    spectrum = full_spectrum(value, height, data_bits)
    bands = local_bands(value, height, data_bits)
    domains, terminal, packing = band_checks(spectrum, bands, height, laminar=False)
    assert len(spectrum) == 5**height
    assert len(domains) == 3**height
    assert rank(spectrum) == 3 * height
    assert packing == len(terminal) == 1 << height
    return {
        "height": height,
        "true_sparsity": len(spectrum),
        "linear_rank": rank(spectrum),
        "distinct_activity_domains": len(domains),
        "activity_packing": str(packing),
    }


def average_basis_margins(ground, target_rank, augmentation):
    bases = [subset for subset in itertools.combinations(ground, target_rank) if rank(value | augmentation for value in subset) == target_rank]
    assert bases
    return {value: Fraction(sum(value in subset for subset in bases), len(bases)) for value in ground}


def unit_cost_checks():
    def cost(margins, marked):
        return sum((value + int(label == marked))**2 for label, value in margins.items()) + int(marked not in margins)

    def boundary_value(point):
        w, v, y, z = [(point >> bit) & 1 for bit in range(4)]
        return 1 - 2 * (w & (v ^ (y & z)))

    true_boundary = walsh(boundary_value(point) for point in range(16))
    line, flat = {0, 1}, {2, 3, 6, 7, 10, 11, 14, 15}
    assert set(true_boundary) == line | flat
    assert sum(map(abs, true_boundary.values())) == 48
    boundary_costs = {}
    for marked in sorted(true_boundary):
        high = average_basis_margins(sorted(flat - {marked}), 4, 16)
        low = {point: Fraction(1, len(line - {marked})) for point in line - {marked}}
        boundary_costs[marked] = str(cost(high | low, marked))
        assert all(
            rank(value | 16 for value in subset + (line_point,)) == 5
            for line_point in low
            for subset in itertools.combinations(high, 4)
            if rank(value | 16 for value in subset) == 4
        )
    assert max(map(Fraction, boundary_costs.values())) == 4
    full_flat_costs = {}
    for dimension in (3, 4):
        number = 1 << dimension
        base = average_basis_margins(tuple(range(number)), dimension + 1, number)
        assert set(base.values()) == {Fraction(dimension + 1, number)}
        calculated = cost(base, 0)
        assert calculated <= 4
        full_flat_costs[dimension] = str(calculated)
    return {"ten_point_deleted_label_costs": boundary_costs, "complete_flat_uniform_costs": full_flat_costs}


def main():
    started = time.monotonic()
    result = {
        "schema": "a28-r074-activation-packing-v1",
        "sharp_trees": [independent_tree(2, False), independent_tree(3, False), independent_tree(2, True)],
        "nested_common_addresses": [nested_common_addresses(height) for height in (2, 3)],
        "crossing_common_addresses": [crossing_common_addresses(height) for height in (2, 3, 4)],
        "unit_cost_checks": unit_cost_checks(),
        "python": sys.version,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "exact integers and fractions; no randomness",
        "elapsed_seconds": time.monotonic() - started,
    }
    Path(sys.argv[1]).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
