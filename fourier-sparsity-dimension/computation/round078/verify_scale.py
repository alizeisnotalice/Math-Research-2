"""Integer audit of the R49 maximum-leaf common-query chain."""
import argparse
import json
import sys
from collections import Counter, defaultdict

LINES = ((0, 1), (2, 4), (3, 7), (5, 6))
TREE_AT_ADDRESS = {point: tree for tree, line in enumerate(LINES) for point in line}


def private_bit(leaves, tree, node):
    return 1 << (3 + tree * (2 * leaves - 1) + node - 1)


def submasks(mask):
    value = 0
    while True:
        yield value
        value = (value - mask) & mask
        if not value:
            break


def parity(mask):
    return mask.bit_count() % 2


def spectral_seed(leaves):
    spectrum = {}
    for tree, (start, end) in enumerate(LINES):
        normals = [alpha for alpha in range(8) if not parity(alpha & (start ^ end))]
        for terminal in range(leaves, 2 * leaves):
            path_mask = right_answers = 0
            node = terminal
            while node > 1:
                parent_bit = private_bit(leaves, tree, node // 2)
                path_mask |= parent_bit
                if node % 2:
                    right_answers |= parent_bit
                node //= 2
            for address_frequency in normals:
                for ancestor_frequency in submasks(path_mask):
                    label = private_bit(leaves, tree, terminal) | ancestor_frequency | address_frequency
                    coefficient = (-1) ** parity((address_frequency & start) | (ancestor_frequency & right_answers))
                    assert label not in spectrum
                    spectrum[label] = coefficient
    return spectrum


def truth_value(leaves, point):
    tree = TREE_AT_ADDRESS[point & 7]
    node = 1
    while node < leaves:
        node = 2 * node + bool(point & private_bit(leaves, tree, node))
    return -1 if point & private_bit(leaves, tree, node) else 1


def walsh(values):
    stride = 1
    while stride < len(values):
        for offset in range(0, len(values), 2 * stride):
            for index in range(offset, offset + stride):
                left, right = values[index], values[index + stride]
                values[index], values[index + stride] = left + right, left - right
        stride *= 2
    return values


def packed_restriction(spectrum, queries, fixed):
    coefficients = defaultdict(int)
    for label, sign in spectrum.items():
        coefficients[label & ~queries] += -sign if parity(label & fixed) else sign
    return {label: value for label, value in coefficients.items() if value}


def projected_multiplicities(spectrum, queries):
    return Counter(label & ~queries for label in spectrum)


def rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length()
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


def leftmost_path(leaves, tree, node):
    mask = 0
    remainder = []
    while node < leaves:
        mask |= private_bit(leaves, tree, node)
        remainder.append((tree, 2 * node + 1))
        node *= 2
    return mask | private_bit(leaves, tree, node), remainder


def size_of(leaves, subtree):
    _, node = subtree
    return leaves >> (node.bit_length() - 1)


def assignment_to_subtree(leaves, subtree):
    tree, node = subtree
    fixed = LINES[tree][0]
    while node > 1:
        if node % 2:
            fixed |= private_bit(leaves, tree, node // 2)
        node //= 2
    return fixed


def exhaustive_leaves(leaves, queries, dimension, numerator, spectrum):
    """For each common leaf, independently read the truth table and run a full FWHT."""
    complement = ((1 << dimension) - 1) ^ queries
    free = list(submasks(complement))
    union = set()
    maximum = 0
    for answer in submasks(queries):
        integers = walsh([truth_value(leaves, answer | x) for x in free])
        actual = {alpha: coeff for alpha, coeff in zip(free, integers) if coeff}
        sparse = packed_restriction(spectrum, queries, answer)
        assert {alpha: coeff * numerator for alpha, coeff in actual.items()} == {
            alpha: coeff * len(free) for alpha, coeff in sparse.items()
        }
        union.update(actual)
        maximum = max(maximum, len(actual))
    return maximum, union


def expected_multiplicity(leaves, active, removed, phase):
    total = 4 * leaves
    if phase == 0:
        return Counter({1: total**2})
    histogram = Counter()
    histogram[removed * total] += 1
    if phase == 1:
        histogram[2] += 6 * leaves**2
    for subtree in active:
        count = size_of(leaves, subtree)
        histogram[total // count] += count**2
    return histogram


def verify_height(height):
    leaves = 1 << (height - 1)
    total = 4 * leaves
    dimension = 2 * total - 1
    spectrum = spectral_seed(leaves)
    assert len(spectrum) == total**2 and rank(spectrum) == dimension
    queries = removed = phase = 0
    records = []
    active = []
    all_terminal_bits = sum(
        private_bit(leaves, tree, terminal)
        for tree in range(4)
        for terminal in range(leaves, 2 * leaves)
    )
    while True:
        quotient = projected_multiplicities(spectrum, queries)
        current_size = len(quotient)
        histogram = Counter(quotient.values())
        assert histogram == expected_multiplicity(leaves, active, removed, phase)
        if phase:
            assert quotient[0] == total * removed
        if current_size <= total:
            break
        if phase == 0:
            mask, active = leftmost_path(leaves, 0, 1)
            expected_m = total**2
            decrease = (29 * leaves**2 - 2) // 3
            batch, fixed = mask | 6, 0
        elif phase == 1:
            mask, remainder = leftmost_path(leaves, 1, 1)
            active.extend([(2, 1), (3, 1), *remainder])
            expected_m = 4 * leaves**2
            decrease = (11 * leaves**2 + 1) // 3
            batch, fixed = mask | 1, 2
        else:
            subtree = min(active, key=lambda item: (-size_of(leaves, item), item))
            count = size_of(leaves, subtree)
            expected_m = count**2
            decrease = (2 * count**2 + 1) // 3
            fixed = assignment_to_subtree(leaves, subtree)
            batch, remainder = leftmost_path(leaves, *subtree)
            active.remove(subtree)
            active.extend(remainder)
        assert not queries & batch
        assert fixed & ~queries == 0
        selected = packed_restriction(spectrum, queries, fixed)
        assert len(selected) == expected_m
        assert all(parity(label & all_terminal_bits) for label in selected)
        certificate_rank = ((expected_m.bit_length() - 1) // 2) + 1
        assert batch.bit_count() == certificate_rank
        assert all(abs(coeff) * (1 << (certificate_rank - 1)) == total for coeff in selected.values())
        assert packed_restriction(spectrum, queries | batch, fixed) == {0: total}
        if height <= 2:
            maximum, union = exhaustive_leaves(leaves, queries, dimension, total, spectrum)
            assert maximum == expected_m and union == set(quotient)
        actual_decrease = current_size - len(projected_multiplicities(spectrum, queries | batch))
        assert decrease == actual_decrease
        records.append({
            "L": current_size, "m": expected_m, "delta_L": actual_decrease,
            "kappa": certificate_rank, "query_rank": batch.bit_count(),
            "multiplicity_histogram": dict(sorted(histogram.items())),
        })
        queries |= batch
        removed += 1
        phase = min(phase + 1, 2)

    remaining_rank = rank(quotient)
    assert remaining_rank + queries.bit_count() == dimension
    if leaves == 1:
        assert current_size == 3 and remaining_rank == 2
    else:
        final_blocks = (leaves + 3) // 3
        assert current_size == 1 + 5 * leaves - 3 * final_blocks
        assert remaining_rank == 4 * leaves - 2 * final_blocks
    if height <= 2:
        _, union = exhaustive_leaves(leaves, queries, dimension, total, spectrum)
        assert union == set(quotient)
    return {
        "H": height, "N": total, "s": len(spectrum), "d": dimension,
        "exhaustive_truth_table_leaves": height <= 2,
        "steps": records, "stop_L": current_size, "cleanup_rank": remaining_rank,
        "stop_multiplicity_histogram": dict(sorted(histogram.items())),
    }


def verify_gate_examples():
    """Full FWHT for the examples with simultaneous growing gate and quadratic parameters."""
    examples = []
    for gates in range(1, 5):
        for pairs in range(1, 5):
            dimension = gates + 2 * pairs
            gate_mask = (1 << gates) - 1
            values = []
            for point in range(1 << dimension):
                polynomial = sum(((point >> (gates + 2 * i)) & 3) == 3 for i in range(pairs)) % 2
                values.append(-1 if polynomial and point & gate_mask == gate_mask else 1)
            integers = walsh(values)
            nonzero_data = integers[1 << gates:]
            assert len(set(map(abs, nonzero_data))) == 1
            assert abs(nonzero_data[0]) == 1 << pairs
            assert rank(label for label, coeff in enumerate(integers) if coeff) == dimension
            examples.append({"gate_codim": gates, "quadratic_pairs": pairs, "vars": dimension})
    return examples


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-height", type=int, default=7)
    args = parser.parse_args()
    print(json.dumps({
        "python": sys.version.split()[0],
        "heights": [verify_height(height) for height in range(1, args.max_height + 1)],
        "gate_examples": verify_gate_examples(),
    }, indent=2))
