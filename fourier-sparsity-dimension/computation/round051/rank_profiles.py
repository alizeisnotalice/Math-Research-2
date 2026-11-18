#!/usr/bin/env python3
"""Exact spectrum intersection profiles and explicit optimal subspaces."""

import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time


LINES = ((0, 1), (2, 4), (3, 7), (5, 6))
ADDRESS_AXES = (0, 1, 6, 7)


def parity(word):
    return word.bit_count() & 1


def span(basis):
    values = [0]
    for row in basis:
        values += [value ^ row for value in values]
    return values


def coordinate_basis(mask):
    return [1 << bit for bit in range(mask.bit_length()) if mask & (1 << bit)]


def combine_spectra(parent, child):
    best = [-1] * (len(parent) + len(child) - 1)
    for left_rank, left_mass in enumerate(parent):
        for right_rank, right_mass in enumerate(child):
            total = left_rank + right_rank
            mass = left_mass + right_mass
            if mass > best[total]:
                best[total] = mass
    return best


def complete_tables(limit):
    result = {1: {"M": [0, 1], "choices": [None, None], "axes": [0, 1]}}
    leaves = 1
    while leaves < limit:
        previous = result[leaves]
        local = previous["M"]
        local_rank = len(local) - 1
        parent_rank = 2 * local_rank + 1
        best = [-1] * (parent_rank + 1)
        for root in (0, 1):
            for first, first_mass in enumerate(local):
                for second, second_mass in enumerate(local):
                    total = root + first + second
                    mass = (first_mass + second_mass) << root
                    if mass > best[total]:
                        best[total] = mass
        greedy = [(0, 0, 0), (0, 1, 0)] + [
            (1, min(rank - 1, local_rank), max(0, rank - 1 - local_rank))
            for rank in range(2, parent_rank + 1)
        ]
        assert best == [(local[left] + local[right]) << root for root, left, right in greedy]
        axes = [
            root | (previous["axes"][left] << 1) | (previous["axes"][right] << (1 + local_rank))
            for root, left, right in greedy
        ]
        assert all(mask.bit_count() == rank for rank, mask in enumerate(axes))
        leaves *= 2
        result[leaves] = {"M": best, "choices": greedy, "axes": axes}
    return result


@cache
def complete_support(leaves):
    if leaves == 1:
        return (1,)
    child = complete_support(leaves // 2)
    left = tuple(word << 1 for word in child)
    right = tuple(word << leaves for word in child)
    points = left + right
    return points + tuple(word | 1 for word in points)


NORMALS = tuple({q for q in range(8) if not parity(q & (left ^ right))} for left, right in LINES)


def normalized_address_weights(address_space):
    return tuple(len(set(address_space) & normals) for normals in NORMALS)


def weighted_allocations(profile, weights):
    best = [0]
    for weight in weights:
        next_best = [-1] * (len(best) + len(profile) - 1)
        for earlier, mass in enumerate(best):
            for rank, count in enumerate(profile):
                total = earlier + rank
                candidate = mass + weight * count
                if candidate > next_best[total]:
                    next_best[total] = candidate
        best = next_best
    greedy = []
    capacity = len(profile) - 1
    order = sorted(range(len(weights)), key=lambda index: (-weights[index], index))
    for total in range(len(best)):
        remaining = total
        assigned = [0] * len(weights)
        for index in order:
            assigned[index] = min(remaining, capacity)
            remaining -= assigned[index]
        assert remaining == 0
        greedy.append(tuple(assigned))
    assert best == [sum(weight * profile[rank] for weight, rank in zip(weights, ranks))
                    for ranks in greedy]
    return best, greedy


def four_line_table(local):
    local_rank = len(local["M"]) - 1
    maximum_rank = 3 + 4 * local_rank
    best = [-1] * (maximum_rank + 1)
    choices = [None] * (maximum_rank + 1)
    for address_rank, axes in enumerate(ADDRESS_AXES):
        weights = normalized_address_weights(span(coordinate_basis(axes)))
        counts, allocations = weighted_allocations(local["M"], weights)
        for private_rank, value in enumerate(counts):
            total = address_rank + private_rank
            if value > best[total]:
                best[total] = value
                choices[total] = (address_rank, *allocations[private_rank])
    axes = []
    for address, *private in choices:
        mask = ADDRESS_AXES[address]
        for group, dimension in enumerate(private):
            mask |= local["axes"][dimension] << (3 + local_rank * group)
        axes.append(mask)
    assert all(mask.bit_count() == rank for rank, mask in enumerate(axes))
    return {"M": best, "choices": choices, "axes": axes}


def four_line_support(leaves):
    block = complete_support(leaves)
    block_rank = 2 * leaves - 1
    return tuple((point << (3 + block_rank * group)) | address
                 for group, normals in enumerate(NORMALS) for point in block for address in normals)


def walsh(values):
    result = list(values)
    width = 1
    while width < len(result):
        for start in range(0, len(result), 2 * width):
            for index in range(start, start + width):
                a, b = result[index], result[index + width]
                result[index], result[index + width] = a + b, a - b
        width *= 2
    return result


def complete_value(leaves, point):
    if leaves == 1:
        return 1 - 2 * (point & 1)
    return complete_value(leaves // 2, point >> (leaves if point & 1 else 1))


def direct_support(leaves, four_lines):
    block_rank = 2 * leaves - 1
    if not four_lines:
        values = (complete_value(leaves, point) for point in range(1 << block_rank))
        denominator = leaves
    else:
        address_to_group = {point: group for group, line in enumerate(LINES) for point in line}
        values = (
            complete_value(leaves, point >> (3 + block_rank * address_to_group[point & 7]))
            for point in range(1 << (3 + 4 * block_rank))
        )
        denominator = 4 * leaves
    transformed = walsh(values)
    assert all(not coefficient or abs(coefficient) * denominator == len(transformed)
               for coefficient in transformed)
    return {frequency for frequency, coefficient in enumerate(transformed) if coefficient}


def subspaces(dimension):
    """RREF enumeration, ordered by least-significant pivot."""
    for rank in range(dimension + 1):
        for pivots in combinations(range(dimension), rank):
            nonpivots = set(range(dimension)) - set(pivots)
            entries = [(row, column) for row, pivot in enumerate(pivots)
                       for column in sorted(nonpivots) if column > pivot]
            for mask in range(1 << len(entries)):
                rows = [1 << column for column in pivots]
                for bit, (row, column) in enumerate(entries):
                    if (mask >> bit) & 1:
                        rows[row] |= 1 << column
                yield rank, rows


def point_mask(points):
    return sum(1 << point for point in points)


def xor_translate(mask, point):
    shifted = 0
    while mask:
        lowbit = mask & -mask
        shifted |= 1 << ((lowbit.bit_length() - 1) ^ point)
        mask ^= lowbit
    return shifted


def small_exhaustion(points, dimension, expected):
    ground = point_mask(points)
    generated = {1}
    for point in points:
        generated |= {mask | xor_translate(mask, point) for mask in generated if not (mask >> point) & 1}
    generated_best = [0] * (dimension + 1)
    generated_counts = Counter()
    for mask in generated:
        rank = mask.bit_count().bit_length() - 1
        generated_best[rank] = max(generated_best[rank], (mask & ground).bit_count())
        generated_counts[rank] += 1
    assert generated_best == expected

    all_best = [0] * (dimension + 1)
    all_counts = Counter()
    seen = set()
    for rank, basis in subspaces(dimension):
        mask = point_mask(span(basis))
        assert mask not in seen
        seen.add(mask)
        all_best[rank] = max(all_best[rank], (mask & ground).bit_count())
        all_counts[rank] += 1
    assert all_best == generated_best
    assert generated <= seen
    return {"point_generated_subspaces_by_rank": dict(sorted(generated_counts.items())),
            "all_subspaces_by_rank": dict(sorted(all_counts.items())),
            "all_subspaces": len(seen), "point_generated_subspaces": len(generated)}


def check_witnesses(support, table, scale):
    counts = table["M"]
    assert len(support) == scale ** 2 and len(set(support)) == len(support)
    for rank, (mass, axes) in enumerate(zip(counts, table["axes"])):
        assert axes.bit_count() == rank
        assert sum(point & axes == point for point in support) == mass
        if rank < len(counts) - 1:
            assert 2 * mass <= scale * rank
    return [rank for rank, mass in enumerate(counts[:-1]) if 2 * mass == scale * rank]


def other_address_check(local, profile):
    ground = [0] * len(profile)
    number = 0
    for rank, basis in subspaces(3):
        weights = normalized_address_weights(span(basis))
        current = [0]
        for weight in weights:
            current = combine_spectra(current, [weight * value for value in local])
        for private, value in enumerate(current):
            ground[rank + private] = max(ground[rank + private], value)
        number += 1
    assert ground == profile and number == 16


def serializable(table):
    return {"M": table["M"], "choices": table["choices"],
            "optimal_coordinate_masks_hex": [hex(value) for value in table["axes"]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--complete-limit", type=int, default=256)
    parser.add_argument("--four-line-limit", type=int, default=64)
    args = parser.parse_args()
    assert all(n >= 4 and n & (n - 1) == 0 for n in (args.complete_limit, args.four_line_limit))
    assert args.four_line_limit <= args.complete_limit
    started = time.perf_counter()
    tables = complete_tables(args.complete_limit)
    complete = {}
    four = {}
    for leaves, table in tables.items():
        points = complete_support(leaves)
        equality = check_witnesses(points, table, leaves)
        complete[leaves] = {"proper_density_equalities": equality, **serializable(table)}
        if leaves <= 4:
            assert set(points) == direct_support(leaves, False)
            complete[leaves]["exhaustion"] = small_exhaustion(points, 2 * leaves - 1, table["M"])
        if leaves >= 4:
            previous = complete[leaves // 2]["proper_density_equalities"]
            assert equality == [0] + [leaves + rank for rank in previous]
        if leaves > args.four_line_limit:
            continue
        four_table = four_line_table(table)
        four_points = four_line_support(leaves)
        outer_equality = check_witnesses(four_points, four_table, 4 * leaves)
        assert outer_equality == [0] + [6 * leaves + rank for rank in equality]
        four[leaves] = {"proper_density_equalities": outer_equality, **serializable(four_table)}
        if leaves <= 4:
            other_address_check(table["M"], four_table["M"])
        if leaves <= 2:
            assert set(four_points) == direct_support(leaves, True)
        if leaves == 1:
            four[leaves]["exhaustion"] = small_exhaustion(four_points, 7, four_table["M"])

    result = {
        "schema": "a28-r051-rank-intersection-profiles-v1",
        "complete_choices_order": ["root_axis_added", "left_rank", "right_rank"],
        "four_line_choices_order": ["address_rank", "first_rank", "second_rank", "third_rank", "fourth_rank"],
        "address_axes_by_rank": [hex(axes) for axes in ADDRESS_AXES],
        "complete": complete,
        "four_lines": four,
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers; no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    print(f"verified every witness through B={args.complete_limit} and four-line B={args.four_line_limit}")


if __name__ == "__main__":
    main()
