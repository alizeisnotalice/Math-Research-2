"""Verify the exact selector-layer identity on deterministic and seeded Boolean trees."""

import json
import platform
import random


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for first in range(0, len(values), 2 * width):
            for position in range(first, first + width):
                a, b = values[position], values[position + width]
                values[position], values[position + width] = a + b, a - b
        width *= 2
    return values


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            index = row.bit_length() - 1
            if index not in pivots:
                pivots[index] = row
                break
            row ^= pivots[index]
    return len(pivots)


def assign_nodes(tree, counter):
    if isinstance(tree, int):
        return tree
    node = counter[0]
    counter[0] += 1
    return node, assign_nodes(tree[0], counter), assign_nodes(tree[1], counter)


def evaluate(tree, selector, data):
    while not isinstance(tree, int):
        node, zero, one = tree
        tree = one if (selector >> node) & 1 else zero
    return 1 - 2 * ((tree >> data) & 1)


def enumerate_leaves(tree):
    if isinstance(tree, int):
        return {tree}
    _, zero, one = tree
    return enumerate_leaves(zero) | enumerate_leaves(one)


def layered_counts(tree, ndata, scale, depth=0):
    if isinstance(tree, int):
        mean = [scale * evaluate(tree, 0, y) for y in range(ndata)]
        return mean, 0, 0, 0, 0
    _, zero, one = tree
    left, lc, lmarks, leff, lzero = layered_counts(zero, ndata, scale, depth + 1)
    right, rc, rmarks, reff, rzero = layered_counts(one, ndata, scale, depth + 1)
    average = [(a + b) // 2 for a, b in zip(left, right)]
    difference = [(a - b) // 2 for a, b in zip(left, right)]
    assert all((a + b) % 2 == 0 for a, b in zip(left, right))
    count = sum(coefficient != 0 for coefficient in walsh(difference))
    effective = len(enumerate_leaves(tree)) > 1
    if effective and (leff == 0 or reff == 0):
        assert count
    return (average, (1 << depth) * count + lc + rc,
            (1 << depth) * int(count > 0) + lmarks + rmarks,
            int(effective) + leff + reff,
            int(effective and count == 0) + lzero + rzero)


def check(tree, data_dimension):
    counter = [0]
    numbered = assign_nodes(tree, counter)
    nodes = counter[0]
    selector_size, data_size = 1 << nodes, 1 << data_dimension
    spectrum = walsh(
        evaluate(numbered, selector, data)
        for data in range(data_size) for selector in range(selector_size)
    )
    support = [frequency for frequency, coefficient in enumerate(spectrum) if coefficient]
    mean, nonzero_layers, marked_cost, effective, unmarked_effective = layered_counts(
        numbered, data_size, selector_size
    )
    mean_support = sum(coefficient != 0 for coefficient in walsh(mean))
    assert len(support) == mean_support + nonzero_layers
    selector_support = {frequency & (selector_size - 1) for frequency in support}
    assert rank(selector_support) == effective
    if effective:
        assert (effective + 1) ** 2 <= 4 * marked_cost
    labels = {frequency >> nodes for frequency in support}
    max_leaf_sparsity = max(
        sum(coefficient != 0 for coefficient in walsh(evaluate(leaf, 0, y)
                                                      for y in range(data_size)))
        for leaf in enumerate_leaves(numbered)
    )
    assert len(labels) ** 2 <= max_leaf_sparsity * len(support)
    return {
        "private_selectors": nodes,
        "effective_selectors": effective,
        "unmarked_effective_selectors": unmarked_effective,
        "marked_cost": marked_cost,
        "true_sparsity": len(support),
        "data_labels": len(labels),
        "max_leaf_sparsity": max_leaf_sparsity,
    }


def main():
    seed = 410038
    rng = random.Random(seed)
    data_dimension = 3
    full_depth_three = lambda leaves: (
        ((leaves[0], leaves[1]), (leaves[2], leaves[3])),
        ((leaves[4], leaves[5]), (leaves[6], leaves[7])),
    )
    constant, opposite = 0, 255
    deterministic = [
        ((constant, opposite), (constant, opposite)),
        (constant, (opposite, constant)),
        ((constant, 3), (252, opposite)),
        (((constant, opposite), (constant, opposite)),
         ((constant, opposite), (constant, opposite))),
    ]
    random_trees = [
        full_depth_three([rng.randrange(256) for _ in range(8)]) for _ in range(32)
    ]
    reports = [check(tree, data_dimension) for tree in deterministic + random_trees]
    result = {
        "python": platform.python_version(),
        "seed": seed,
        "fixed_data_dimension": data_dimension,
        "all_checks_passed": len(reports),
        "deterministic_examples": reports[:len(deterministic)],
        "random_example_maximum_sparsity": max(row["true_sparsity"] for row in reports),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
