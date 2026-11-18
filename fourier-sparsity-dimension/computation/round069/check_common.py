"""Integer Walsh checks for dependent shared labels and near-Sidon private bags."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform


def sign(bits):
    return 1 - 2 * (bits.bit_count() & 1)


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


def walsh(values):
    width = 1
    while width < len(values):
        for offset in range(0, len(values), 2 * width):
            for pos in range(offset, offset + width):
                left, right = values[pos], values[pos + width]
                values[pos], values[pos + width] = left + right, left - right
        width *= 2
    return values


def circuit(powers, shift):
    return [1 << (shift + j) for j in range(powers)] + [((1 << powers) - 1) << shift]


def check_code(lengths):
    offset = 0
    common = []
    relations = []
    for length in lengths:
        current = circuit(length - 1, offset)
        assert rank(current) == length - 1
        assert all(sign(col) == -1 for col in current)
        mask = ((1 << length) - 1) << len(common)
        relations.append(mask)
        offset += length - 1
        common.extend(current)
    common_rank = rank(common)
    assert common_rank == offset
    assert len(common) - common_rank == len(relations) == rank(relations)
    assert len({a ^ b for a, b in combinations(common, 2)}) == len(common) * (len(common) - 1) // 2
    return common, common_rank


def check_matchings(common, rounds):
    size = len(common)
    odd = size - 1
    all_pairs = set()
    differences = set()
    for row in range(rounds):
        pairs = [(odd, row)] + [
            ((row + delta) % odd, (row - delta) % odd)
            for delta in range(1, size // 2)
        ]
        assert {point for pair in pairs for point in pair} == set(range(size))
        for i, j in pairs:
            pair = tuple(sorted((i, j)))
            difference = common[i] ^ common[j]
            assert pair not in all_pairs
            assert difference not in differences
            all_pairs.add(pair)
            differences.add(difference)
    # The two-label common sector of each actual signed Boolean plane is checked.
    plane = walsh([1, 1, 1, -1])
    assert plane == [2, 2, 2, -2]
    return len(all_pairs), len(differences)


def shared_levels(case):
    if case == "four_common":
        common_rank = 4
        labels = [1 << j for j in range(common_rank)]
        levels = [sum(sign(x & alpha) for alpha in labels) // 2 for x in range(1 << common_rank)]
        expected_histogram = {0: 6, 1: 8, 2: 2}
    else:
        common_rank = 9
        labels = circuit(common_rank, 0)
        levels = [
            (sum(sign(x & alpha) for alpha in labels[:-1]) - sign(x & labels[-1])) // 4
            for x in range(1 << common_rank)
        ]
        expected_histogram = {0: 252, 1: 240, 2: 20}
    histogram = Counter(map(abs, levels))
    assert dict(histogram) == expected_histogram
    return common_rank, labels, levels, histogram


def verify_local(case):
    common_rank, labels, levels, counts = shared_levels(case)
    size_common = 1 << common_rank
    values = []
    for z in range(32):
        for level in levels:
            if abs(level) == 2:
                value = level // 2
            elif abs(level) == 1:
                value = level * (-1 if (z & 3) == 0 else 1)
            else:
                value = sign(((z >> 4) & 1) ^ (((z >> 2) & 1) * ((z >> 3) & 1)))
            values.append(value)

    spectrum = {index: coeff for index, coeff in enumerate(walsh(values[:])) if coeff}
    groups = Counter(index >> common_rank for index in spectrum)
    assert {index for index in spectrum if index < size_common} == set(labels)
    expected_support = {0, 1, 2, 3, 16, 20, 24, 28}
    assert set(groups) == expected_support
    multiplicity = Counter(a ^ b for a, b in combinations(groups, 2))
    for a, b, c in combinations(groups, 3):
        d = a ^ b ^ c
        if d in groups:
            assert all(multiplicity[diff] == 2 for diff in (a ^ b, a ^ c, b ^ c))

    for x, level in enumerate(levels):
        actual = {index: coeff for index, coeff in enumerate(walsh(values[x::size_common]))
                  if coeff}
        target = (set([0]) if abs(level) == 2 else set(range(4)) if abs(level) == 1
                  else set([16, 20, 24, 28]))
        assert set(actual) == target
        assert actual.get(0, 0) == 16 * level

    if case == "four_common":
        assert len(spectrum) == 60
        assert set(groups.values()) == {4, 8}
    cost = len(labels) + (
        Fraction(4, 3) * counts[1] + Fraction(9, 4) * counts[0]
    ) / size_common
    assert cost == (Fraction(529, 96) if case == "four_common" else Fraction(6007, 512))
    return {
        "local_case": case,
        "common_rank": common_rank,
        "simultaneous_common_labels": len(labels),
        "private_projection_labels": len(groups) - 1,
        "private_fibers": dict(sorted(groups.items())),
        "local_spectrum_size": len(spectrum),
        "absolute_mean_counts_scaled_by_2": dict(sorted(counts.items())),
        "local_budget": str(cost),
        "truth_table_entries": len(values),
    }


def verify_families(n, local_four, local_ten):
    lengths = [12, 20] * n
    common, common_rank = check_code(lengths)
    length = len(common)
    match_rounds = length // 2
    ports, differences = check_matchings(common, match_rounds)
    assert common_rank == 30 * n
    four_bags = length // 4

    ten_labels, ten_rank = check_code([10] * n)
    assert ten_rank == 9 * n
    return {
        "n": n,
        "distance_12": {
            "common_labels": length,
            "common_rank": common_rank,
            "dependency_rank": len(lengths),
            "exact_dependency_distance": min(lengths),
            "matching_subfamily": {
                "bags": match_rounds,
                "ports": ports,
                "distinct_common_directions": differences,
            },
            "four_common_subfamily": {
                "bags": four_bags,
                "total_data_rank": common_rank + 5 * four_bags,
                "sparsity": local_four["local_spectrum_size"] * four_bags ** 2,
            },
        },
        "distance_10": {
            "common_labels": len(ten_labels),
            "common_rank": ten_rank,
            "dependency_rank": n,
            "exact_dependency_distance": 10,
            "bags": n,
            "total_data_rank": ten_rank + 5 * n,
            "sparsity": local_ten["local_spectrum_size"] * n ** 2,
        },
    }


if __name__ == "__main__":
    four = verify_local("four_common")
    ten = verify_local("ten_common_odd_signs")
    print(json.dumps({
        "python": platform.python_version(),
        "local_cases": [four, ten],
        "family_checks": [verify_families(n, four, ten) for n in (1, 2, 4)],
    }, indent=2))
