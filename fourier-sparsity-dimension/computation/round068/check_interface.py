"""Exact sparse verification of two four-point blocks with a large shared span."""
import itertools
import json
import platform


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


def sign(value):
    return 1 - 2 * (value.bit_count() & 1)


def verify(a):
    m = 1 << a
    d = m - 1
    normal_dim = a + d
    columns = [
        (1 << j) | (1 << ((j + 1) % d)) | (1 << ((j + 3) % d))
        for j in range(d)
    ]
    assert rank(columns) == d
    simplex = [
        sum((sign(eta & point) == -1) << (point - 1) for point in range(1, m))
        for eta in range(m)
    ]

    low_images = [0] + columns + [left ^ right for left, right in itertools.combinations(columns, 2)]
    expected_nonzero_differences = m * len(low_images) - 1
    min_transformed_weight = min(
        (codeword ^ vector).bit_count()
        for eta, codeword in enumerate(simplex)
        for index, vector in enumerate(low_images)
        if eta or index
    )
    assert min_transformed_weight >= 3
    assert all(codeword.bit_count() == m // 2 for codeword in simplex[1:])

    first = set()
    second = set()
    # These are the unnormalized integer numerators of the original address spectrum.
    signed_address = {}
    for eta in range(m):
        for cell in range(m):
            data = 0 if cell == 0 else 1 << (cell - 1)
            p_data = 0 if cell == 0 else columns[cell - 1]
            alpha = eta | (data << a)
            first.add(alpha | (1 << normal_dim))
            second.add(eta | ((p_data ^ simplex[eta]) << a) | (2 << normal_dim))
            signed_address[alpha] = sign(eta & cell)
    assert len(first) == len(second) == m * m
    assert not first.intersection(second)
    assert rank(first) == rank(second) == normal_dim + 1
    global_rank = rank(first | second)
    assert global_rank == normal_dim + 2

    # The inverse transform in the address variables recovers the one-hot cells.
    for point in range(m):
        for cell in range(m):
            data = 0 if cell == 0 else 1 << (cell - 1)
            assert sum(sign(eta & point) * signed_address[eta | (data << a)] for eta in range(m)) == (
                m if point == cell else 0
            )

    if a == 5:
        # Check all pair differences, and verify Booleanity directly by sparse convolution.
        left_diff = {left ^ right for left, right in itertools.combinations(first, 2)}
        right_diff = {left ^ right for left, right in itertools.combinations(second, 2)}
        assert len(left_diff) == len(right_diff) == expected_nonzero_differences
        assert left_diff.isdisjoint(right_diff)
        square = {0: m * m}
        for (alpha, coeff), (beta, other) in itertools.combinations(signed_address.items(), 2):
            delta = alpha ^ beta
            square[delta] = square.get(delta, 0) + 2 * coeff * other
        assert {frequency: count for frequency, count in square.items() if count} == {0: m * m}


    guards = [1 << d, 2 << d]
    private = [4 << d, 8 << d]
    data_cells = [0] + [1 << j for j in range(d)]
    character_blocks = [
        {guard ^ atom ^ (axis if c else 0) for atom in atoms for c in (0, 1)}
        for guard, axis, atoms in zip(guards, private, (data_cells, [0] + columns))
    ]
    assert all(len(block) == 2 * m for block in character_blocks)
    character_global = rank(set.union(*character_blocks))
    character_interface = sum(rank(block) for block in character_blocks) - character_global
    assert character_global == m + 3 and character_interface == m - 1
    difference_sets = [
        {left ^ right for left, right in itertools.combinations(block, 2)}
        for block in character_blocks
    ]
    assert difference_sets[0].isdisjoint(difference_sets[1])
    core_size = 4 * m
    singleton_coefficients = [sign(core_point & core_frequency)
                              for core_point in range(core_size)
                              for core_frequency in range(core_size)]
    assert all(abs(coefficient) == 1 for coefficient in singleton_coefficients)

    return {
        "a": a,
        "cells": m,
        "spectral_points_per_block": m * m,
        "difference_points_per_block": expected_nonzero_differences,
        "min_transformed_data_difference_weight": min_transformed_weight,
        "interface_rank": 2 * (normal_dim + 1) - global_rank,
        "fixed_slice_rank": normal_dim + 1,
        "global_projected_rank": global_rank,
        "relation_rank": 2 * (m * m - m - a),
        "original_tail_sparsity": 4 * m * m,
        "dense_core_tail_sparsity": 8 * m * m,
        "full_pair_difference_and_sparse_square_checks": a == 5,
        "rank_one_slice_variant": {
            "labels_per_block": 2 * m,
            "interface_rank": character_interface,
            "projected_rank": character_global,
            "max_slice_rank": 1,
            "tail_sparsity": len(singleton_coefficients),
        },
    }


if __name__ == "__main__":
    print(json.dumps({"python": platform.python_version(), "checks": [verify(5), verify(7)]}, indent=2))
