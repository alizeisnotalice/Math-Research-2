"""Exact finite checks for arbitrary-frontier amplification."""
import itertools
import json
import random
import sys
from fractions import Fraction


def optimal_cost(count):
    base = 1 << (count.bit_length() - 1)
    return base * (3 * count - 2 * base)


def support(values):
    values = values[:]
    stride = 1
    while stride < len(values):
        for offset in range(0, len(values), 2 * stride):
            for index in range(offset, offset + stride):
                left, right = values[index], values[index + stride]
                values[index], values[index + stride] = left + right, left - right
        stride *= 2
    return [index for index, coefficient in enumerate(values) if coefficient]


def rank(vectors):
    basis = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = vector
                break
            vector ^= basis[pivot]
    return len(basis)


def check_binary_trees(max_leaves):
    costs = [0, 1]
    for count in range(2, max_leaves + 1):
        costs.append(2 * min(costs[left] + costs[count - left] for left in range(1, count)))
        assert costs[-1] == optimal_cost(count)


def address_partitions(remaining):
    if not remaining:
        yield []
        return
    block = remaining
    least_point = remaining & -remaining
    while block:
        if block & least_point:
            for later in address_partitions(remaining ^ block):
                yield [block, *later]
        block = (block - 1) & remaining


def check_partitions():
    supports = [support([(mask >> point) & 1 for point in range(8)]) for mask in range(256)]
    partition_count = 0
    optimum_count = 0
    for blocks in address_partitions(255):
        count = len(blocks)
        tally = sum(len(supports[mask]) for mask in blocks)
        assert tally >= optimal_cost(count)
        partition_count += 1
        if tally == optimal_cost(count):
            assert all(mask.bit_count() * len(supports[mask]) == 8 for mask in blocks)
            assert rank(itertools.chain.from_iterable(supports[mask] for mask in blocks)) <= count - 1
            optimum_count += 1
    return partition_count, optimum_count


def lookup_tables(counts):
    leaves = [()]
    frontiers = {1: leaves[:]}
    for width in range(2, counts[-1] + 1):
        leaf = min(leaves, key=lambda path: (len(path), path))
        leaves.remove(leaf)
        leaves.extend((leaf + (0,), leaf + (1,)))
        frontiers[width] = leaves[:]

    nodes = sorted({leaf[:depth] for leaf in leaves for depth in range(len(leaf))})
    node_index = {node: index for index, node in enumerate(nodes)}
    data_index = {}
    for layer, count in enumerate(counts):
        for cell in frontiers[count]:
            data_index[layer, cell] = len(nodes) + len(data_index)

    lookups = []
    for address in range(1 << len(nodes)):
        leaf = ()
        while leaf in node_index:
            leaf += ((address >> node_index[leaf]) & 1,)
        columns = []
        for layer, count in enumerate(counts):
            cell = next(path for path in frontiers[count] if leaf[:len(path)] == path)
            columns.append(data_index[layer, cell])
        lookups.append(columns)
    return lookups


def check_nested_walsh(seed):
    rng = random.Random(seed)
    tests = 0
    for counts in itertools.combinations_with_replacement(range(1, 5), 3):
        lookups = lookup_tables(counts)
        dimension = counts[-1] - 1 + sum(counts)
        for _ in range(2):
            while True:
                values = [rng.choice((-1, 1)) for _ in range(8)]
                old_spectrum = support(values)
                sizes = [
                    sum(vector.bit_length() == layer for vector in old_spectrum)
                    for layer in range(1, 4)
                ]
                if rank(old_spectrum) == 3 and min(sizes) >= 1:
                    break
            lifted_values = []
            for vector in range(1 << dimension):
                columns = lookups[vector & (len(lookups) - 1)]
                point = sum(((vector >> index) & 1) << layer for layer, index in enumerate(columns))
                lifted_values.append(values[point])
            spectrum = support(lifted_values)
            expected_size = int(0 in old_spectrum) + sum(
                layer_size * optimal_cost(count) for layer_size, count in zip(sizes, counts)
            )
            assert len(spectrum) == expected_size
            assert rank(spectrum) == dimension
            tests += 1
    return tests


def check_rounding():
    tests = 0
    for counts in itertools.combinations_with_replacement(range(1, 17), 4):
        bases = [1 << (count.bit_length() - 1) for count in counts]
        probabilities = [Fraction(count, base) - 1 for count, base in zip(counts, bases)]
        endpoints = sorted({Fraction(0), Fraction(1), *probabilities})
        first_moments = [0] * len(counts)
        second_moments = [0] * len(counts)
        for left, right in zip(endpoints, endpoints[1:]):
            trial = [
                (2 if (left + right) / 2 < probability else 1) * base
                for probability, base in zip(probabilities, bases)
            ]
            assert trial == sorted(trial)
            for layer, rounded_count in enumerate(trial):
                first_moments[layer] += (right - left) * rounded_count
                second_moments[layer] += (right - left) * rounded_count ** 2
        assert first_moments == list(counts)
        assert second_moments == list(map(optimal_cost, counts))
        tests += 1
    return tests


if __name__ == "__main__":
    check_binary_trees(512)
    partition_count, optimum_count = check_partitions()
    seed = 74
    summary = {
        "python": sys.version.split()[0],
        "seed": seed,
        "binary_tree_cost_B_range": [1, 512],
        "F2_3_partitions": partition_count,
        "cost_minimizers": optimum_count,
        "nested_walsh_tests": check_nested_walsh(seed),
        "monotone_rounding_four_tuples": check_rounding(),
    }
    print(json.dumps(summary, indent=2))
