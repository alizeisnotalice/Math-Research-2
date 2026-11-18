#!/usr/bin/env python3
"""Independent integer audit of the five-dimensional second-weight proof."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
from time import perf_counter


SHAPES = {
    "point": (0,),
    "line": (0, 1),
    "three": (0, 1, 2),
    "plane": (0, 1, 2, 3),
    "three_times_two": (0, 1, 2, 3, 4, 5),
    "sidon": (0, 1, 2, 4, 8, 15),
}
MASS_LIMITS = (Fraction(1, 2), Fraction(7, 8), Fraction(5, 4), Fraction(2))


def span(vectors):
    result = {0}
    for vector in vectors:
        result |= {point ^ vector for point in result}
    return frozenset(result)


def subspaces(dimension):
    return {space for basis in combinations(range(1, 16), dimension)
            if len(space := span(basis)) == 1 << dimension}


def character_masks(active):
    return tuple(sum(1 << index for index, point in enumerate(active)
                     if (point & frequency).bit_count() % 2)
                 for frequency in range(16))


def abs_walsh(active, mask):
    return tuple(abs(len(active) - 2 * (mask ^ character).bit_count())
                 for character in character_masks(active))


def zero_set_census(zeros, spaces):
    active = tuple(index for index in range(16) if index not in zeros)
    n = len(active)
    masks = character_masks(active)
    counts = Counter()
    for row in range(1 << (n - 1)):
        spectrum = tuple(abs(n - 2 * (row ^ mask).bit_count()) for mask in masks)
        if max(spectrum[1:]) >= 8 or sum(spectrum[1:]) >= 32:
            continue
        if any(sum(spectrum[index] for index in labels) >= 14 for labels in spaces[2]):
            continue
        if any(sum(spectrum[index] for index in labels) >= 20 for labels in spaces[3]):
            continue
        counts[spectrum[0]] += 1
    return dict(sorted(counts.items()))


def transform_label(point):
    # Invertible mixing x0 <- x0+x1, x2 <- x2+x3 followed by a translation.
    return point ^ ((point >> 1) & 1) ^ ((point >> 1) & 4) ^ 10


def check_sidon_partitions(spaces):
    sidon = set(SHAPES["sidon"])
    assert {a ^ b for a, b in combinations(sidon, 2)} == set(range(1, 16))
    allowed_sizes = {0: {1}, 1: {2}, 2: {2, 3}, 3: {2}, 4: {6}}
    profiles = {}
    for dimension, subspaces_by_dimension in spaces.items():
        profiles[dimension] = Counter()
        for direction in subspaces_by_dimension:
            classes = {frozenset(sidon & {point ^ h for h in direction}) for point in sidon}
            profile = tuple(sorted(map(len, classes)))
            assert set(profile) & allowed_sizes[dimension]
            if dimension == 3:
                assert profile == (2, 4)
            profiles[dimension][",".join(map(str, profile))] += 1
    return profiles


def check_projection_comparison(census, spaces):
    bounds = {}
    slope_ranks = {
        "point": (4,),
        "line": (3, 4),
        "three": (4,),
        "plane": (3,),
        "three_times_two": (3,),
        "sidon": (4,),
    }
    for name, zeros in SHAPES.items():
        absolute_coefficients = abs_walsh(zeros, 0)
        weights = tuple(Fraction(value, len(zeros)) for value in absolute_coefficients)
        upper = Fraction(0)
        for value in sorted(set(weights[1:]), reverse=True):
            if not value:
                continue
            labels = {label for label in range(1, 16) if weights[label] >= value}
            direction = frozenset(labels | {0})
            dimension = (len(direction) - 1).bit_length()
            assert direction in spaces[dimension]
            lower_values = [weight for weight in weights[1:] if weight < value]
            upper += (value - max(lower_values, default=0)) * MASS_LIMITS[dimension - 1]
        for rank in slope_ranks[name]:
            lower = (1 << (rank - 2)) * Fraction(min(census[name]), 16) ** 2
            assert lower > upper
            bounds[f"{rank}:{name}"] = {"lower_over_q": str(lower), "upper_over_q": str(upper)}
    return bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = perf_counter()
    spaces = {dimension: subspaces(dimension) for dimension in range(5)}
    assert [len(spaces[dimension]) for dimension in range(5)] == [1, 15, 35, 15, 1]
    census_spaces = {dimension: tuple(tuple(space - {0}) for space in members)
                     for dimension, members in spaces.items()}
    census = {name: zero_set_census(zeros, census_spaces) for name, zeros in SHAPES.items()}
    transformed = {
        name: zero_set_census(tuple(map(transform_label, zeros)), census_spaces)
        for name, zeros in SHAPES.items()
    }
    assert census == transformed
    reference_bytes = args.reference.read_bytes()
    reference = json.loads(reference_bytes)
    assert reference["three_half"] == {name: {str(key): count for key, count in result.items()}
                                       for name, result in census.items()}
    projection = check_projection_comparison(census, spaces)
    result = {
        "schema": "a28-r038-independent-second-weight-audit-v1",
        "parameter_census": census,
        "rows_per_complete_census": sum(1 << (15 - len(zeros)) for zeros in SHAPES.values()),
        "affine_transform_recheck": transformed == census,
        "core_partitions": check_sidon_partitions(spaces),
        "normalized_projection_comparisons": projection,
        "reference_sha256": sha256(reference_bytes).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers and rational numbers; no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"checked {2 * result['rows_per_complete_census']} parameter rows and all core partitions")


if __name__ == "__main__":
    main()
