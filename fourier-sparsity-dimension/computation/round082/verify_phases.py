"""Exact Walsh and binary-pencil checks for multibit overlap potentials."""
import itertools
import json
import sys
from collections import Counter


def walsh(values):
    coefficients = list(values)
    stride = 1
    while stride < len(coefficients):
        for offset in range(0, len(coefficients), 2 * stride):
            for index in range(offset, offset + stride):
                left, right = coefficients[index], coefficients[index + stride]
                coefficients[index], coefficients[index + stride] = left + right, left - right
        stride *= 2
    return coefficients


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


def affstats(values):
    size = len(values)
    support = [i for i, coefficient in enumerate(walsh(values)) if coefficient]
    return len(support), rank(label | size for label in support)


def alternating_rows(size, mask):
    rows = [0] * size
    for index, (left, right) in enumerate(itertools.combinations(range(size), 2)):
        if (mask >> index) & 1:
            rows[left] ^= 1 << right
            rows[right] ^= 1 << left
    return rows


def check_rank_two_pencils():
    dimension = 4
    rows = [alternating_rows(dimension, i) for i in range(1 << 6)]
    ranks = list(map(rank, rows))
    eligible = tuple(i for i, r in enumerate(ranks) if r <= 2)
    count = Counter()
    sharp = Counter()
    for m0, m1, m2 in itertools.product(eligible, repeat=3):
        m3 = m0 ^ m1 ^ m2
        if ranks[m3] > 2:
            continue
        matrices = [m0, m1, m2, m3]
        address_rank = rank(vector for i in matrices for vector in rows[i])
        sparsity = sum(1 << ranks[i] for i in matrices)
        assert 4 * (4 + address_rank) ** 2 <= 16 * sparsity
        distinct = len(set(matrices))
        contains_zero = str(0 in matrices)
        count[(distinct, contains_zero)] += 1
        if 4 * (4 + address_rank) ** 2 == 16 * sparsity:
            sharp[(distinct, contains_zero)] += 1
    return (
        [{"distinct_matrices": d, "zero_present": z == "True", "count": n}
         for (d, z), n in sorted(count.items())],
        [{"distinct_matrices": d, "zero_present": z == "True", "count": n}
         for (d, z), n in sorted(sharp.items())],
    )


def single_flip(size):
    return tuple(-1 if x == 0 else 1 for x in range(1 << size))


def equality_example(data_dimension):
    address_dimension = data_dimension + 2
    seed = single_flip(data_dimension)
    values = []
    for point in range(1 << (data_dimension + address_dimension)):
        addresses = point & ((1 << address_dimension) - 1)
        data = point >> address_dimension
        active = addresses & 1
        shift = (addresses >> 2) if active else 0
        sign = (-1) ** (active and (addresses & 2) != 0)
        values.append(sign * seed[data ^ shift])
    return tuple(values)


def check_equality():
    checked = []
    for data in range(2, 7):
        values = equality_example(data)
        sparsity, dimension = affstats(values)
        assert (sparsity, dimension) == (4 * (1 << data), 2 * data + 3)
        slices = [affstats(values[address::1 << (data + 2)]) for address in range(1 << (data + 2))]
        assert set(slices) == {(1 << data, data + 1)}
        assert (dimension + 1) ** 2 * (1 << data) == (data + 2) ** 2 * sparsity
        checked.append({
            "data_dimension": data, "R_parent": dimension,
            "s_parent": sparsity, "all_address_slices": len(slices),
        })
    return checked


def check_variable_support():
    # Data and address bits in each private group, plus one group-specific guard.
    data = 2
    group_address_width = data + 2
    innermost = equality_example(data)
    group_bits = 1 + 2 * data + 2
    low_mask = (1 << group_bits) - 1
    values = []
    for x in range(1 << (1 + 2 * group_bits)):
        group = x & 1
        group_value = (x >> (1 + group_bits * group)) & low_mask
        sign = (-1) ** (group_value & 1)
        values.append(sign * innermost[group_value >> 1])
    sparsity, dimension = affstats(values)
    assert (sparsity, dimension) == (64, 15)
    assert (dimension + 1) ** 2 == 4 * sparsity
    seed = single_flip(data)
    repetitions = Counter()
    for z0, u_left, u_right in itertools.product(range(2), range(16), range(16)):
        # Put the two private guards in the data cube together with two data blocks.
        group_u = (u_left, u_right)[z0]
        branch = []
        for x in range(1 << (2 * (data + 1))):
            private = (x >> (z0 * (data + 1))) & ((1 << (data + 1)) - 1)
            sign = (-1) ** (private & 1)
            branch.append(sign * innermost[group_u + ((private >> 1) << group_address_width)])
        assert affstats(branch) == (4, 3)
        for beta, coefficient in enumerate(walsh(branch)):
            if coefficient:
                repetitions[beta] += 1
                assert abs(coefficient) * 2 == len(branch)
    assert len(repetitions) == 8 and set(repetitions.values()) == {256}
    return {"s_parent": sparsity, "R_parent": dimension, "address_branches": 512,
            "labels": len(repetitions), "repetitions_per_label": min(repetitions.values())}


def check_counterexamples():
    results = []
    line_count = 0
    for dimension in range(2, 8):
        base = single_flip(dimension)
        seed = tuple(
            (-1) ** ((x >> 1) & 1) if x & 1 else base[x >> 2]
            for x in range(1 << (2 + dimension))
        )
        parent = tuple(
            seed[(x >> 2) ^ (2 if x & 3 == 3 else 0)]
            for x in range(1 << (4 + dimension))
        )
        child_s, child_r = affstats(seed)
        parent_s, parent_r = affstats(parent)
        assert (child_s, child_r) == ((1 << (dimension + 1)) + 2, dimension + 3)
        assert (parent_s, parent_r) == ((1 << (dimension + 1)) + 8, dimension + 5)
        assert (parent_r + 1) ** 2 * child_s > (child_r + 1) ** 2 * parent_s
        for addresses in itertools.combinations(range(4), 2):
            values = tuple(
                parent[((x >> 1) << 2) | addresses[x & 1]]
                for x in range(1 << (3 + dimension))
            )
            assert affstats(values) == (child_s, child_r)
            line_count += 1
        assert {
            affstats(parent[u::4]) for u in range(4)
        } == {(child_s, child_r)}
        results.append({"dimension": dimension, "parent_s": parent_s, "parent_R": parent_r,
                        "child_s": child_s, "child_R": child_r})
    return results, line_count


if __name__ == "__main__":
    pencil_counts, sharp_counts = check_rank_two_pencils()
    counterex, affine_lines = check_counterexamples()
    print(json.dumps({
        "python": sys.version.split()[0],
        "pencil_counts": pencil_counts,
        "sharp_pencils": sharp_counts,
        "full_support_equality": check_equality(),
        "variable_support_equality": check_variable_support(),
        "nonfull_support_counterexamples": counterex,
        "counterexample_single_bit_affine_restrictions": affine_lines,
    }, indent=2))
