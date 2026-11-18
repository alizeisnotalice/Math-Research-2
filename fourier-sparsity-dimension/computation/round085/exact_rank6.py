"""Finite integer certificate for four rank-two images spanning six dimensions."""
from collections import Counter
import itertools
import json
import sys

SIZE = 6
PAIRS = tuple(itertools.combinations(range(SIZE), 2))


def points(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def wedge(x, y):
    result = 0
    for index, (left, right) in enumerate(PAIRS):
        if (((x >> left) & (y >> right)) ^ ((x >> right) & (y >> left))) & 1:
            result |= 1 << index
    return result


def rows(mask):
    result = [0] * SIZE
    for index in points(mask):
        left, right = PAIRS[index]
        result[left] ^= 1 << right
        result[right] ^= 1 << left
    return result


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


def extend(subspace, vectors):
    for vector in vectors:
        if not ((subspace >> vector) & 1):
            subspace |= sum(1 << (vector ^ x) for x in points(subspace))
    return subspace


def subspaces(dimension):
    full = (1 << (1 << dimension)) - 1
    current = {1}
    for dim in range(dimension + 1):
        following = set()
        for space in sorted(current):
            yield space, dim
            outside = full ^ space
            members = tuple(points(space))
            while outside:
                x = (outside & -outside).bit_length() - 1
                coset = sum(1 << (x ^ y) for y in members)
                following.add(space | coset)
                outside ^= coset
        current = following


def profiles():
    matrices = {}
    for x in range(1, 1 << SIZE):
        for y in range(x + 1, 1 << SIZE):
            matrices.setdefault(wedge(x, y), (x, y))
    assert len(matrices) == 651
    ranks = tuple(rank(rows(mask)) for mask in range(1 << len(PAIRS)))
    first, second = wedge(1, 2), wedge(4, 8)
    initial = (1 << 16) - 1
    answers = {}
    total = 0
    for third, third_vectors in sorted(matrices.items()):
        extended = extend(initial, third_vectors)
        for fourth, fourth_vectors in sorted(matrices.items()):
            if extend(extended, fourth_vectors).bit_count() != 64:
                continue
            masks = [0]
            for matrix in (first, second, third, fourth):
                masks.extend(old ^ matrix for old in masks[:])
            profile = tuple(ranks[mask] for mask in masks)
            answers.setdefault(profile, (first, second, third, fourth))
            total += 1
    assert total == 354176
    assert len(answers) == 62
    return answers, total


def cut_certificate(profile, spaces):
    zeros = sum(1 << i for i, value in enumerate(profile) if value == 0)
    rank_two = sum(1 << i for i, value in enumerate(profile) if value == 2)
    candidates = [(space, outside & rank_two, payoff)
                  for space, outside, payoff in spaces if not outside & zeros]
    candidates.sort(key=lambda item: -item[2])
    for size in (2, 3):
        for collection in itertools.combinations(candidates, size):
            if sum(item[2] for item in collection) < 6:
                continue
            lower_rank_covered = 0
            for _, selected, _ in collection:
                if lower_rank_covered & selected:
                    break
                lower_rank_covered |= selected
            else:
                return [space for space, _, _ in collection]
    raise AssertionError(profile)


def verify_cost(profile, certificate, space_dimensions):
    total_payoff = sum(5 - space_dimensions[space] for space in certificate)
    assert total_payoff >= 6 and len(certificate) <= 3
    for value, polar_rank in enumerate(profile):
        tests = sum(not ((space >> value) & 1) for space in certificate)
        assert tests <= (1 << (polar_rank // 2)) - 1
        # This is 2^(r/2)(1−2^(−r)) >= (5/4)*tests, multiplied by 4·2^r.
        assert 4 * (1 << (polar_rank // 2)) * ((1 << polar_rank) - 1) >= (
            5 * tests * (1 << polar_rank)
        )


if __name__ == "__main__":
    examples, total = profiles()
    space_dimensions = dict(subspaces(4))
    full = (1 << 16) - 1
    spaces = [(space, full ^ space, 5 - dim)
              for space, dim in sorted(space_dimensions.items()) if dim < 4]
    records = []
    for profile, example in sorted(examples.items()):
        certificate = cut_certificate(profile, spaces)
        verify_cost(profile, certificate, space_dimensions)
        records.append(dict(profile=profile, representative_forms=example, cut_subspaces=certificate))
    counts = Counter(len(row["cut_subspaces"]) for row in records)
    print(json.dumps({"python":sys.version.split()[0], "ordered_plane_pairs":total,
                      "profile_count":len(records), "cut_count_histogram":counts,
                      "certificates":records}, indent=2))
