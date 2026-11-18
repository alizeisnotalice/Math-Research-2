#!/usr/bin/env python3
"""Exact cover contractions and Boolean graph spectra for A28 round 59."""

import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
import sys
import time


def parity(value):
    return value.bit_count() & 1


def bits(mask):
    while mask:
        least = mask & -mask
        yield least.bit_length() - 1
        mask ^= least


def rank_binary(vectors):
    basis = {}
    for value in vectors:
        while value and value.bit_length() in basis:
            value ^= basis[value.bit_length()]
        if value:
            basis[value.bit_length()] = value
    return len(basis)


def rank_real_integer(columns):
    basis = {}
    for column in columns:
        value = list(column)
        for index, coordinate in enumerate(value):
            if not coordinate:
                continue
            prior = basis.get(index)
            if prior is None:
                basis[index] = value
                break
            value[index:] = [
                prior[index] * current - coordinate * old
                for current, old in zip(value[index:], prior[index:])
            ]
            gcd = math.gcd(*value[index:])
            if gcd > 1:
                value[index:] = [entry // gcd for entry in value[index:]]
    return len(basis)


def linear_subspaces(n):
    subspaces = {frozenset({0})}
    for vector in range(1, 1 << n):
        for subspace in tuple(subspaces):
            if vector not in subspace:
                subspaces.add(subspace | {vector ^ entry for entry in subspace})
    return sorted(subspaces, key=lambda subspace: (len(subspace), sorted(subspace)))


def affine_masks(n, subspaces):
    return {
        sum(1 << (point ^ translate) for point in subspace)
        for subspace in subspaces
        for translate in range(1 << n)
    }


def affine_data(n, mask):
    points = tuple(bits(mask))
    representative = points[0]
    normal = tuple(
        frequency
        for frequency in range(1, 1 << n)
        if all(not parity(frequency & (point ^ representative)) for point in points)
    )
    assert len(points) * (1 + len(normal)) == 1 << n
    column = {
        frequency: 1 - 2 * parity(frequency & representative)
        for frequency in normal
    }
    return normal, column


def positive_integer_covers(n, target, masks):
    size = 1 << n
    choices = [
        sorted(mask for mask in masks if next(bits(mask)) == point)
        for point in range(size)
    ]

    def recurse(point, remaining, selected):
        if point == size:
            assert not any(remaining)
            yield selected
            return
        for batch in itertools.combinations_with_replacement(choices[point], remaining[point]):
            updated = list(remaining)
            for mask in batch:
                for covered in bits(mask):
                    updated[covered] -= 1
            if min(updated) >= 0:
                assert updated[point] == 0
                yield from recurse(point + 1, updated, selected + batch)

    yield from recurse(0, [target] * size, ())


def check_cover_contractions(n, targets):
    subspaces = linear_subspaces(n)
    masks = affine_masks(n, subspaces)
    domains = {mask: affine_data(n, mask) for mask in masks}
    profiles = []
    for target in targets:
        count = overlaps = strict = equalities = 0
        for cover in positive_integer_covers(n, target, masks):
            count += 1
            distinct = tuple(sorted(set(cover)))
            overlaps += sum(mask.bit_count() for mask in distinct) > 1 << n
            normals = tuple(freq for mask in distinct for freq in domains[mask][0])
            for subspace in subspaces:
                outside = [freq for freq in range(1, 1 << n) if freq not in subspace]
                rank_before = rank_binary(subspace)
                quotient_rank = rank_binary((*subspace, *normals)) - rank_before
                real_rank = rank_real_integer(
                    [domains[mask][1].get(freq, 0) for freq in outside]
                    for mask in distinct
                )
                assert quotient_rank <= real_rank <= len(distinct) - 1
                strict += quotient_rank < real_rank
                equalities += quotient_rank == real_rank
        profiles.append({
            "dimension": n,
            "target_weight": target,
            "cover_count": count,
            "overlapping_covers": overlaps,
            "quotient_subspaces": len(subspaces),
            "strict_inequalities": strict,
            "equalities": equalities,
        })
    return profiles


def field_products(n, polynomial):
    size = 1 << n

    def multiply(left, right):
        result = 0
        for index in range(n):
            if (right >> index) & 1:
                result ^= left
            left <<= 1
            if left & size:
                left ^= polynomial
        return result

    products = [[multiply(x, y) for y in range(size)] for x in range(size)]
    for x in range(1, size):
        assert len(set(products[x])) == size
    return products


def walsh(values):
    width = 1
    while width < len(values):
        for start in range(0, len(values), width * 2):
            for left in range(start, start + width):
                right = left + width
                first, second = values[left], values[right]
                values[left], values[right] = first + second, first - second
        width *= 2
    return {frequency: value for frequency, value in enumerate(values) if value}


def check_field_graph(n, polynomial):
    size = 1 << n
    products = field_products(n, polynomial)
    z_size = size // 2
    pairings = [
        sum((products[lam][1 << index] & 1) << index for index in range(n))
        for lam in range(size)
    ]
    assert len(set(pairings)) == size
    values = [
        1 - 2 * (y == products[u][z])
        for z in range(z_size)
        for y in range(size)
        for u in range(size)
    ]
    spectrum = walsh(values)
    scale = len(values)
    expected = {0: scale - (scale >> (n - 1))}
    activity_masks = set()
    frequency_normals = set()
    unscaled_columns = []
    for lam in range(1, size):
        activities = [0] * z_size
        for u in range(size):
            y_freq = pairings[products[lam][u]]
            omega = y_freq & (z_size - 1)
            activities[omega] |= 1 << u
        assert all(activity.bit_count() == 2 for activity in activities)
        for omega, mask in enumerate(activities):
            assert mask not in activity_masks
            activity_masks.add(mask)
            normal, column = affine_data(n, mask)
            frequency_normals.update(normal)
            unscaled_columns.append(column)
            data_label = pairings[lam] | (omega << n)
            for frequency, sign in {0: 1, **column}.items():
                expected[frequency | (data_label << n)] = -(1 << (n + 1)) * sign
    assert len(activity_masks) == math.comb(size, 2)
    assert spectrum == expected
    assert len(spectrum) == 1 + z_size**2 * (size - 1)
    assert rank_binary(spectrum) == 3 * n - 1
    assert rank_binary(frequency_normals) == n
    conditional_ranks = []
    for old_rank in range(n + 1):
        filter_threshold = 1 << old_rank
        real_rank = rank_real_integer(
            [column.get(frequency, 0) for frequency in range(filter_threshold, size)]
            for column in unscaled_columns
        )
        assert real_rank == size - filter_threshold
        assert n - old_rank <= real_rank
        conditional_ranks.append(real_rank)
    return {
        "n": n,
        "polynomial_binary": bin(polynomial),
        "truth_table_size": scale,
        "distinct_activity_domains": 1 + len(activity_masks),
        "conditional_real_ranks_on_coordinate_flag": conditional_ranks,
        "normal_rank": n,
        "sparsity": len(spectrum),
        "global_rank": rank_binary(spectrum),
    }


def main():
    started = time.monotonic()
    covers = check_cover_contractions(2, (1, 2, 3)) + check_cover_contractions(3, (1,))
    graphs = [check_field_graph(n, polynomial) for n, polynomial in enumerate((7, 11, 19, 37, 67), 2)]
    result = {
        "schema": "a28-r059-affine-cover-rank-v1",
        "covers": covers,
        "boolean_field_graphs": graphs,
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers; no random inputs",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.monotonic() - started,
    }
    Path(sys.argv[1]).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
