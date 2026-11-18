#!/usr/bin/env python3
"""Exact minimum parity decision tree census for all four-bit truth tables."""

import argparse
from collections import Counter, defaultdict
from functools import cache
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time


def query_masks(dimension):
    size = 1 << dimension
    return tuple(sum(((point & parity).bit_count() & 1) << point for point in range(size))
                 for parity in range(size))


def query_partitions(dimension):
    size = 1 << dimension
    partitions = []
    for parity in range(1, size):
        pivot_bit = 1 << (parity.bit_length() - 1)
        lower_bits = pivot_bit - 1
        points_without_pivot = [(point & lower_bits) | ((point & ~lower_bits) << 1)
                                for point in range(size // 2)]
        partitions.append(tuple(
            tuple(point | (pivot_bit * (((point & parity).bit_count() & 1) ^ branch))
                  for point in points_without_pivot)
            for branch in (0, 1)
        ))
    return tuple(partitions)


def section(mask, points):
    return sum(((mask >> point) & 1) << index for index, point in enumerate(points))


def opt_table(dimension, previous):
    leaf = [0] * (1 << (1 << dimension))
    partitions = query_partitions(dimension)
    full = len(leaf) - 1
    for mask in range(1, full):
        leaf[mask] = min(1 + previous[section(mask, left)] + previous[section(mask, right)]
                         for left, right in partitions)
    return leaf


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


def spectrum(mask, characters, size):
    return [size - 2 * (mask ^ character).bit_count() for character in characters]


def spectral_profile(mask, characters, size):
    active = [query for query, value in enumerate(spectrum(mask, characters, size)) if value]
    origin = active[0]
    return len(active), rank(active), rank(frequency ^ origin for frequency in active)


def affine_ambient_check(mask, dimension):
    """Independent state recurrence keeping physical point sets, not compressed rows."""
    characters = query_masks(dimension)
    full = (1 << (1 << dimension)) - 1

    @cache
    def optimum(domain):
        if domain & mask in (0, domain):
            return 0
        return min(costs(domain).values())

    def costs(domain):
        result = {}
        for query, character in enumerate(characters[1:], 1):
            one = domain & character
            zero = domain ^ one
            if zero and one:
                result[query] = 1 + optimum(zero) + optimum(one)
        return result

    def tree(domain):
        if domain & mask == 0:
            return {"sign": 1}
        if domain & mask == domain:
            return {"sign": -1}
        query = next(query for query, cost in costs(domain).items() if cost == optimum(domain))
        one = domain & characters[query]
        return {"query_mask": query, "zero": tree(domain ^ one), "one": tree(one)}

    def evaluate(node, point):
        if "sign" in node:
            return node["sign"]
        branch = "one" if (point & node["query_mask"]).bit_count() & 1 else "zero"
        return evaluate(node[branch], point)

    def node_count(node):
        if "sign" in node:
            return 0
        return 1 + node_count(node["zero"]) + node_count(node["one"])

    witness = tree(full)
    expected = [1 - 2 * ((mask >> point) & 1) for point in range(1 << dimension)]
    assert [evaluate(witness, point) for point in range(1 << dimension)] == expected
    assert node_count(witness) == optimum(full)
    return optimum(full), costs(full), witness


def dimension_profile(dimension, table):
    characters = query_masks(dimension)
    size = 1 << dimension
    counts = defaultdict(Counter)
    representatives = {}
    maximal_tree_size = max(table)
    for mask, best in enumerate(table):
        profile = spectral_profile(mask, characters, size)
        counts[best][profile] += 1
        if best == maximal_tree_size:
            representatives.setdefault(profile, mask)
    return counts, representatives


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check-mask", type=lambda value: int(value, 0), action="append", default=[])
    args = parser.parse_args()
    started = time.perf_counter()
    dimension = 4
    tables = [[0, 0]]
    for index in range(1, dimension + 1):
        tables.append(opt_table(index, tables[-1]))
    table = tables[dimension]
    counts, representatives = dimension_profile(dimension, table)
    maximum = max(table)

    witnesses = []
    mask_set = set(representatives.values()) | set(args.check_mask)
    partitions = query_partitions(dimension)
    for mask in sorted(mask_set):
        assert 0 <= mask < len(table)
        chars = query_masks(dimension)
        check, root_costs, witness = affine_ambient_check(mask, dimension)
        assert check == table[mask]
        compressed_costs = {
            query: 1 + tables[dimension - 1][section(mask, left)]
            + tables[dimension - 1][section(mask, right)]
            for query, (left, right) in enumerate(partitions, 1)
        }
        assert root_costs == compressed_costs
        s, linear, affine = spectral_profile(mask, chars, 1 << dimension)
        witnesses.append({
            "truth_mask_hex": hex(mask),
            "negative_points": [point for point in range(1 << dimension) if mask >> point & 1],
            "minimum_internal_nodes": check,
            "minimum_total_nodes": 2 * check + 1,
            "sparsity": s,
            "linear_rank": linear,
            "affine_rank": affine,
            "unnormalized_walsh": spectrum(mask, chars, 1 << dimension),
            "fixed_root_internal_nodes_by_query": root_costs,
            "optimal_tree": witness,
        })
    output = {
        "schema": "a28-r047-parity-decision-tree-census-v1",
        "dimension": dimension,
        "truth_tables_by_dimension": [len(row) for row in tables],
        "histogram_by_dimension": [dict(sorted(Counter(row).items())) for row in tables],
        "joint_profile_order": ["sparsity", "linear_rank", "affine_rank"],
        "dimension_four_joint_counts": {
            size: {",".join(map(str, profile)): count for profile, count in sorted(profiles.items())}
            for size, profiles in sorted(counts.items())
        },
        "hardest_internal_nodes": maximum,
        "hardest_count": sum(counts[maximum].values()),
        "hardest_counts_by_joint_profile": {
            ",".join(map(str, profile)): count for profile, count in sorted(counts[maximum].items())
        },
        "witnesses": witnesses,
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers; no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    assert sum(map(sum, (profile.values() for profile in counts.values()))) == len(table)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"checked {len(table)} four-bit functions: max {maximum}, attained {output['hardest_count']} times")


if __name__ == "__main__":
    main()
