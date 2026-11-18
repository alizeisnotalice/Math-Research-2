"""Exact Walsh and GF(2) checks for reverse-order and nonparallel gluing."""
from collections import Counter
from itertools import combinations
import json
import platform


def rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            index = vector.bit_length() - 1
            if index not in pivots:
                pivots[index] = vector
                break
            vector ^= pivots[index]
    return len(pivots)


def augmented_rank(points):
    return rank((point << 1) | 1 for point in points)


def restrict(spectrum, query, answer):
    pivot = 1 << (query.bit_length() - 1)
    result = {}
    for frequency, coefficient in spectrum.items():
        if frequency & pivot:
            frequency ^= query
            if answer:
                coefficient = -coefficient
        result[frequency] = result.get(frequency, 0) + coefficient
    return {frequency: coefficient for frequency, coefficient in result.items() if coefficient}


def compare_one_step(parent, left, right, common_rank):
    """Sign of h+2sqrt(s0) compared with 4sqrt(s1); this test has s1=s2."""
    assert len(left) == len(right)
    tail = 16 * len(left) - common_rank ** 2 - 4 * len(parent)
    if tail < 0:
        return 1
    test = 16 * common_rank ** 2 * len(parent) - tail ** 2
    return (test > 0) - (test < 0)


def walsh(values):
    width = 1
    while width < len(values):
        for offset in range(0, len(values), 2 * width):
            for i in range(offset, offset + width):
                a, b = values[i], values[i + width]
                values[i], values[i + width] = a + b, a - b
        width *= 2
    return values


def binary_spectrum(leaves):
    signs = {}
    for leaf in range(leaves, 2 * leaves):
        level = [(1 << (leaf - 1), 1)]
        node = leaf
        while node > 1:
            parent = node // 2
            factor = -1 if node & 1 else 1
            bit = 1 << (parent - 1)
            level += [(frequency ^ bit, factor * coefficient) for frequency, coefficient in level]
            node = parent
        signs.update(level)
    assert len(signs) == leaves ** 2
    return signs


def compare_truth(leaves, spectrum):
    values = []
    size = 1 << (2 * leaves - 1)
    for point in range(size):
        node = 1
        while node < leaves:
            node = 2 * node + ((point >> (node - 1)) & 1)
        values.append(1 - 2 * ((point >> (node - 1)) & 1))
    w = walsh(values)
    assert {i: c for i, c in enumerate(w) if c} == {
        i: c * (size // leaves) for i, c in spectrum.items()
    }


def binary_reverse(leaves, exhaustive):
    original = binary_spectrum(leaves)
    if leaves <= 8:
        compare_truth(leaves, original)
    active = [original]
    total_subtraction = 0
    profile = []
    for j, node in enumerate(range(leaves // 2, leaves)):
        following = []
        j_s = leaves ** 2 - 3 * leaves * j // 2
        j_dim = 2 * leaves - 2 * j
        expected_common_rank = j_dim - 4
        query = 1 << (node - 1)
        for parent in active:
            assert (len(parent), augmented_rank(parent) + 1) == (j_s, j_dim)
            left, right = (restrict(parent, query, answer) for answer in (0, 1))
            common = left.keys() & right.keys()
            assert len(common) == j_s - 2 * leaves
            assert all(left[x] == right[x] for x in common)
            intersection_rank = augmented_rank(left) + augmented_rank(right) - augmented_rank(left.keys() | right.keys())
            assert intersection_rank == augmented_rank(common) == expected_common_rank
            payment_sign = compare_one_step(parent, left, right, intersection_rank)
            assert payment_sign == (-1 if j < leaves // 2 - 1 else 1)
            following.extend((left, right) if exhaustive else (left,))
        active = following
        total_subtraction += (1 << j) * expected_common_rank
        profile.append({
            "fixed_bottom_selectors": j, "sparsity": j_s, "D": j_dim,
            "common_rank": expected_common_rank, "common_support": j_s - 2 * leaves,
            "one_step_sign": payment_sign,
        })
    assert all(len(leaf) == leaves ** 2 // 4 and augmented_rank(leaf) + 1 == leaves for leaf in active)
    assert (1 << (leaves // 2)) * leaves - total_subtraction == 2 * leaves
    return {
        "terminal_characters": leaves, "full_truth_checked": leaves <= 8,
        "exhaustive_resolution": exhaustive,
        "consecutive_failing_levels": leaves // 2 - 1,
        "theoretical_terminal_restrictions": 1 << (leaves // 2),
        "last_checked_leaf_count": len(active),
        "weighted_common_rank_subtraction": total_subtraction,
        "level_profile": profile,
    }


def four_line_tree(depth):
    lines = [(0, 1), (2, 4), (3, 7), (5, 6)]
    next_index = 0
    terminal_mask = 0

    def descend(step):
        nonlocal next_index, terminal_mask
        if step == 0:
            position = 1 << next_index
            next_index += 1
            terminal_mask |= position
            return {position: 1}
        start = next_index
        next_index += 3
        result = {}
        for left, right in lines:
            addresses = [(direction << start, 1 - 2 * ((direction & left).bit_count() & 1))
                         for direction in range(8) if (direction & (left ^ right)).bit_count() % 2 == 0]
            for frequency, coefficient in descend(step - 1).items():
                for address, phase in addresses:
                    result[frequency ^ address] = coefficient * phase
        return result

    coefficients = descend(depth)
    return coefficients, terminal_mask, next_index


def four_line_audit(depth):
    spectrum, terminal_mask, dim = four_line_tree(depth)
    terminal_count = 4 ** depth
    assert dim == rank(spectrum) == 2 * terminal_count - 1
    assert len(spectrum) == terminal_count ** 2
    assert augmented_rank(spectrum) == dim
    pair_counts = Counter(left ^ right for left, right in combinations(spectrum, 2))
    direction_types = Counter()
    for query, collisions in pair_counts.items():
        left, right = [restrict(spectrum, query, answer) for answer in (0, 1)]
        assert len(left) == len(right) == len(spectrum) - 3 * collisions // 2
        common_rank = augmented_rank(left) + augmented_rank(right) - augmented_rank(left.keys() | right.keys())
        assert compare_one_step(spectrum, left, right, common_rank) == -1
        direction_types[(query & terminal_mask).bit_count()] += 1

    assert sum(direction_types.values()) == len(pair_counts)
    query = 1
    left, right = [restrict(spectrum, query, answer) for answer in (0, 1)]
    root_share = augmented_rank(left.keys() & right.keys())
    assert root_share == terminal_count
    return {
        "macro_depth": depth, "terminal_characters": terminal_count, "full_spectrum": len(spectrum),
        "strict_directions_checked": len(pair_counts),
        "strict_directions_by_terminal_coordinate_count": dict(sorted(direction_types.items())),
        "root_scheduled_common_rank": root_share,
    }


if __name__ == "__main__":
    print(json.dumps({
        "python": platform.python_version(),
        "reverse": [binary_reverse(b, b <= 16) for b in (4, 8, 16, 32)],
        "four_lines": [four_line_audit(depth) for depth in (1, 2)],
    }, indent=2))
