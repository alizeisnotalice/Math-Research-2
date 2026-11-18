"""Integer verification of the arbitrary-rank affine-pencil bounds."""
from collections import Counter
from functools import cache
import itertools
import json
import random
import sys

SEED = 8302


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


def walsh(values):
    values = list(values)
    stride = 1
    while stride < len(values):
        for start in range(0, len(values), 2 * stride):
            for index in range(start, start + stride):
                left, right = values[index], values[index + stride]
                values[index], values[index + stride] = left + right, left - right
        stride *= 2
    return values


@cache
def pairs(size):
    return tuple(itertools.combinations(range(size), 2))


def rows(size, mask):
    result = [0] * size
    for index, (first, second) in enumerate(pairs(size)):
        if (mask >> index) & 1:
            result[first] ^= 1 << second
            result[second] ^= 1 << first
    return tuple(result)


def wedge_mask(size, first, second):
    return sum(1 << i for i, pair in enumerate(pairs(size)) if pair == tuple(sorted((first, second))))


def pencil_masks(generators):
    values = [generators[0]]
    for difference in generators[1:]:
        values.extend(value ^ difference for value in values[:])
    return values


def check_pencil(size, generators, test_all_directions):
    v = len(generators) - 1
    count = 1 << v
    masks = pencil_masks(generators)
    matrices = [rows(size, mask) for mask in masks]
    ranks = [rank(matrix) for matrix in matrices]
    a = rank(vector for matrix in matrices for vector in matrix)
    widths = sum(1 << r for r in ranks)
    shift = v + 2
    lhs = count * (shift + a) ** 2
    rhs = shift ** 2 * widths
    assert lhs <= rhs
    if lhs == rhs:
        assert max(ranks) == 0 or (min(ranks) == max(ranks) == 2 and a == shift)
    if max(ranks) >= 4:
        assert 112 * lhs <= 81 * rhs
    difference_ranks = {direction: rank(rows(size, masks[0] ^ masks[direction])) for direction in range(1, count)}
    largest = max(difference_ranks.values())
    selected = difference_ranks if test_all_directions else {
        next(direction for direction, r in difference_ranks.items() if r == largest): largest
    }
    pairs_checked = 0
    for direction, difference_rank in selected.items():
        if difference_rank == 0:
            continue
        assert 2 * count * a <= 2 * count * difference_rank + v * (
            2 * sum(ranks) - count * difference_rank
        )
        for i in range(count):
            j = i ^ direction
            if i >= j:
                continue
            pair_rank = rank(matrices[i] + matrices[j])
            assert 2 * (pair_rank - difference_rank) <= ranks[i] + ranks[j] - difference_rank
            pairs_checked += 1
    mode = "constant" if largest == 0 else (
        "d_ge_mean" if count * largest >= sum(ranks) else "d_lt_mean"
    )
    return max(ranks), mode, pairs_checked


def random_pencils():
    rng = random.Random(SEED)
    result = Counter()
    checked_pairs = 0
    for size in (3, 4, 5, 6, 8, 10, 12, 16):
        for v in range(2, 8):
            for sample in range(16):
                generators = [rng.getrandbits(len(pairs(size))) for _ in range(v + 1)]
                if sample == 0:
                    # Include high-rank constant pencils explicitly.
                    generators[1:] = [0] * v
                maximum, mode, pair_count = check_pencil(
                    size, generators, test_all_directions=v <= 5
                )
                result[maximum, mode] += 1
                checked_pairs += pair_count
    return result, checked_pairs


def eval_quad(mask, address, size):
    parity = 0
    for index, (first, second) in enumerate(pairs(size)):
        if (mask >> index) & 1:
            parity ^= ((address >> first) & 1) & ((address >> second) & 1)
    return parity


def true_parent(size, generators):
    v = len(generators) - 1
    cube = 1 << (size + v)
    signs = []
    for point in range(cube):
        address = point & ((1 << size) - 1)
        x = point >> size
        phase = eval_quad(generators[0], address, size)
        translation = sum(
            eval_quad(mask, address, size) << j
            for j, mask in enumerate(generators[1:])
        )
        source = -1 if x ^ translation == 0 else 1
        signs.append(((-1) ** phase) * source)
    support = [index for index, coeff in enumerate(walsh(signs)) if coeff]
    return len(support), rank(index | cube for index in support)


def true_examples():
    size = 4
    a = wedge_mask(size, 0, 1)
    b = wedge_mask(size, 2, 3)
    c = a ^ wedge_mask(size, 0, 2)
    no_common_line = true_parent(size, (a, a ^ b, a ^ c))
    assert no_common_line == (16, 7)

    size = 5
    a = wedge_mask(size, 0, 1)
    b = wedge_mask(size, 2, 3)
    c = wedge_mask(size, 0, 4)
    gap_sharp = true_parent(size, (a, a ^ b, a ^ c))
    assert gap_sharp == (28, 8)
    assert 112 * (gap_sharp[1] + 1) ** 2 * 4 == 81 * 16 * gap_sharp[0]

    near = []
    for v in range(2, 8):
        size = v + 1
        generators = [0] + [wedge_mask(size, 0, i) for i in range(1, v + 1)]
        sparsity, homogeneous_rank = true_parent(size, generators)
        assert sparsity == 4 * (1 << v) - 3
        assert homogeneous_rank == 2 * v + 2
        near.append(dict(v=v, s=sparsity, R=homogeneous_rank))
    return {
        "all_rank_two_no_common_line_s_R": no_common_line,
        "sharp_gap_s_R": gap_sharp, "one_zero_near_sharp_s_R": near,
    }


if __name__ == "__main__":
    cases, paired_cases = random_pencils()
    print(json.dumps({
        "python": sys.version.split()[0], "seed": SEED,
        "affine_pencils": sum(cases.values()), "paired_rank_inequalities": paired_cases,
        "rank_mode_counts": [
            dict(max_rank=maximum, mode=mode, count=count)
            for (maximum, mode), count in sorted(cases.items())
        ],
        "true_examples": true_examples(),
    }, indent=2))
