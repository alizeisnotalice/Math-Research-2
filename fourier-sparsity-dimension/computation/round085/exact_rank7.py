"""Finite integer certificate for five rank-two images spanning seven dimensions."""
from collections import Counter
from functools import cache
import itertools
import json
import sys

from exact_rank6 import points, rank, subspaces


def dot(x, y):
    return (x & y).bit_count() & 1


def image(columns, value):
    result = 0
    for index in points(value):
        result ^= columns[index]
    return result


@cache
def basis(space):
    vectors = []
    for vector in points(space):
        if rank([*vectors, vector]) > len(vectors):
            vectors.append(vector)
    return tuple(vectors)


def transform(space, columns):
    return sum(1 << image(columns, vector) for vector in points(space))


def coordinate_kernels():
    # Each coordinate projection K -> F_2^2 has codimension at most two.
    kernels = tuple(space for space, dim in subspaces(3) if dim >= 1)
    assert len(kernels) == 15
    index = {space: i for i, space in enumerate(kernels)}
    group = [
        tuple(index[transform(space, columns)] for space in kernels)
        for columns in itertools.permutations(range(1, 8), 3) if rank(columns) == 3
    ]
    assert len(group) == 168
    pair_forms = []
    for space in kernels:
        if space.bit_count() != 2:
            pair_forms.append((0, 0, 0))
            continue
        annihilator = [z for z in range(1, 8) if all(dot(z, x) == 0 for x in points(space))]
        a, b = annihilator[:2]
        pair_forms.append(tuple(
            (b if (a >> coordinate) & 1 else 0) ^ (a if (b >> coordinate) & 1 else 0)
            for coordinate in range(3)
        ))
    return kernels, group, pair_forms


def rank_profiles():
    kernels, group, pair_forms = coordinate_kernels()
    represented = {}
    count = 0
    for combination in itertools.combinations_with_replacement(range(len(kernels)), 5):
        intersections = [(1 << 8) - 1]
        for index in combination:
            intersections.extend(space & kernels[index] for space in intersections[:])
        if any(intersections[31 ^ (1 << i)] != 1 for i in range(5)):
            continue
        if any(tuple(sorted(g[index] for index in combination)) < combination for g in group):
            continue
        forms = [(0, 0, 0)]
        for index in combination:
            forms.extend(tuple(x ^ y for x, y in zip(form, pair_forms[index])) for form in forms[:])
        values = []
        for active, form in enumerate(forms):
            vectors = basis(intersections[31 ^ active])
            square = tuple(
                sum(dot(image(form, x), y) << i for i, y in enumerate(vectors)) for x in vectors
            )
            values.append(2 * active.bit_count() - 2 * len(vectors) + rank(square))
        assert values[0] == 0 and all(0 <= r <= 6 and r % 2 == 0 for r in values)
        assert all(values[1 << i] == 2 for i in range(5))
        represented.setdefault(tuple(values), tuple(kernels[index] for index in combination))
        count += 1
    assert count == 109 and len(represented) == 97
    return represented, count


def search(profile, spaces):
    zeros = sum(1 << i for i, value in enumerate(profile) if value == 0)
    rank_two = sum(1 << i for i, value in enumerate(profile) if value == 2)
    candidates = [(space, outside & rank_two, payoff)
                  for space, outside, payoff in spaces if not outside & zeros]
    candidates.sort(key=lambda item: -item[2])
    for size in (2, 3):
        for collection in itertools.combinations(candidates, size):
            if sum(item[2] for item in collection) < 7:
                continue
            previous = 0
            for _, selected, _ in collection:
                if previous & selected:
                    break
                previous |= selected
            else:
                return tuple(space for space, _, _ in collection)
    raise AssertionError(profile)


def verify(profile, certificate, space_dimensions):
    assert sum(6 - space_dimensions[space] for space in certificate) >= 7
    for label, polar_rank in enumerate(profile):
        count = sum(not ((space >> label) & 1) for space in certificate)
        scale = 1 << (polar_rank // 2)
        assert count <= scale - 1
        assert 4 * scale * ((1 << polar_rank) - 1) >= 5 * count * (1 << polar_rank)


if __name__ == "__main__":
    profiles, count = rank_profiles()
    dimensions = dict(subspaces(5))
    full = (1 << 32) - 1
    spaces = [(space, full ^ space, 6 - dim)
              for space, dim in sorted(dimensions.items()) if dim < 5]
    records = []
    for profile, kernels in sorted(profiles.items()):
        certificate = search(profile, spaces)
        verify(profile, certificate, dimensions)
        records.append(dict(profile=profile, coordinate_kernel_masks=kernels, cut_subspaces=certificate))
    print(json.dumps({"python":sys.version.split()[0], "kernel_orbits":count,
                      "profile_count":len(profiles),
                      "cut_count_histogram":Counter(len(row["cut_subspaces"]) for row in records),
                      "certificates":records}, indent=2))
