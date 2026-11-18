#!/usr/bin/env python3
"""Exact checks of Boolean and general polynomial tail-restriction witnesses."""

import hashlib
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


def linear_subspaces(n):
    subspaces = {frozenset({0})}
    for vector in range(1, 1 << n):
        for subspace in tuple(subspaces):
            if vector not in subspace:
                subspaces.add(subspace | {vector ^ point for point in subspace})
    return subspaces


def coordinate_order(subspace):
    coordinates = [0]
    span = {0}
    for vector in sorted(subspace):
        if vector not in span:
            shifted = [vector ^ point for point in coordinates]
            coordinates.extend(shifted)
            span.update(shifted)
    assert span == subspace
    return coordinates


def affine_flats(n):
    for subspace in linear_subspaces(n):
        coordinates = coordinate_order(subspace)
        unused = set(range(1 << n))
        while unused:
            representative = min(unused)
            points = [representative ^ vector for vector in coordinates]
            unused.difference_update(points)
            yield points


def boolean_value(query, data):
    address = query & 3
    if address:
        return 1 - 2 * ((data >> (address - 1)) & 1)
    return 1 - 2 * (query == 0)


def boolean_witness(a):
    n = a + 2
    spectrum = walsh(
        boolean_value(query, data)
        for data in range(8)
        for query in range(1 << n)
    )
    actual_head = {freq: value for freq, value in spectrum.items() if freq >> n == 0}
    actual_tail = {freq: value for freq, value in spectrum.items() if freq >> n != 0}
    assert set(actual_head) == set(range(1 << n))
    assert set(actual_tail) == {freq | (data << n) for data in (1, 2, 4) for freq in range(4)}
    assert all(value % 32 == 16 for value in actual_head.values())
    assert all(value % 32 == 0 for value in actual_tail.values())
    flats_count = retaining_count = minimum_count = 0
    minimum = len(spectrum) + 1
    counts_by_dimension = {}
    for flat in affine_flats(n):
        flats_count += 1
        tail_profiles = [walsh((point & 3) == address for point in flat) for address in (1, 2, 3)]
        if not all(tail_profiles):
            continue
        retaining_count += 1
        assert {point & 3 for point in flat} == {0, 1, 2, 3}
        assert all(len(profile) == 4 for profile in tail_profiles)
        head = walsh(boolean_value(point, 0) if point & 3 == 0 else 0 for point in flat)
        assert len(head) >= 4
        total = len(head) + sum(map(len, tail_profiles))
        dimension = len(flat).bit_length() - 1
        counts_by_dimension[dimension] = counts_by_dimension.get(dimension, 0) + 1
        if total < minimum:
            minimum = total
            minimum_count = 0
        minimum_count += total == minimum
    assert minimum == 16
    return {
        "a": a,
        "highest_odd_flat_size": len(actual_head),
        "tail_size": len(actual_tail),
        "outer_rank": 3,
        "query_affine_flats": flats_count,
        "rank_preserving_flats": retaining_count,
        "rank_preserving_flats_by_dim": counts_by_dimension,
        "minimum_restricted_sparsity": minimum,
        "minimizing_flats": minimum_count,
    }


def multiplication(a, polynomial):
    size = 1 << a

    def mul(x, y):
        product = 0
        for index in range(a):
            if (y >> index) & 1:
                product ^= x
            x <<= 1
            if x & size:
                x ^= polynomial
        return product

    for x in range(1, size):
        assert len({mul(x, y) for y in range(size)}) == size
    return mul


def spread_witness(a, polynomial):
    size = 1 << a
    mul = multiplication(a, polynomial)
    augmented_blocks = [tuple(x | (mul(slope, x) << a) for x in range(size)) for slope in range(size)]
    augmented_blocks.extend((0, y << a) for y in range(2, size))
    vectors = [vector for block in augmented_blocks for vector in block if vector]
    assert len(vectors) == len(set(vectors)) == size**2 - 2
    assert set(vectors) == set(range(1, size**2)) - {size}
    query_size = size**2 // 2
    forbidden = set()
    total_tail = 0
    for block in augmented_blocks:
        equations = {
            (vector & (size - 1)) | ((vector >> (a + 1)) << a): (vector >> a) & 1
            for vector in block
        }
        assert len(equations) == len(block)
        values = [
            all((frequency & point).bit_count() % 2 == sign for frequency, sign in equations.items())
            for point in range(query_size)
        ]
        spectrum = walsh(values)
        assert set(spectrum) == set(equations)
        assert all(spectrum[frequency] == (1 - 2 * sign) * sum(values) for frequency, sign in equations.items())
        for frequency, side in equations.items():
            if frequency:
                forbidden.add((frequency, 1 - side))
        total_tail += len(spectrum)
    assert len(forbidden) == 2 * (query_size - 1)
    outer_rank = len(augmented_blocks)
    assert total_tail == size**2 + 2 * size - 4
    assert query_size == (total_tail - outer_rank) // 2 + 1
    return {
        "a": a,
        "field_polynomial": bin(polynomial),
        "outer_rank": outer_rank,
        "blocked_affine_hyperplanes": len(forbidden),
        "head_size": query_size,
        "tail_size": total_tail,
        "minimum_restricted_sparsity": total_tail + query_size,
    }


def main():
    started = time.monotonic()
    result = {
        "schema": "a28-r063-tail-restrictions-v1",
        "boolean_witnesses": [boolean_witness(a) for a in (2, 3, 4)],
        "general_polynomial_witnesses": [spread_witness(a, p) for a, p in enumerate((7, 11, 19, 37), 2)],
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers; no randomness",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.monotonic() - started,
    }
    Path(sys.argv[1]).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
