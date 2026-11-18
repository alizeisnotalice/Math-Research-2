#!/usr/bin/env python3
"""Exact finite audit of the R34 phase normal form and its R32 reference."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform
import sys
import time


def parity(value):
    return value.bit_count() & 1


def coordinates(point):
    return tuple((point >> index) & 1 for index in range(7))


def get_branches():
    branches = {"A": [], "B": [], "C": []}
    for point in range(128):
        a, b, c, d, e, f, g = coordinates(point)
        if (a, b, e, f) == (1, 0, 1, 1):
            branches["A"].append(point)
        if (a, b, c, g) == (0, 1, 1, 1):
            branches["B"].append(point)
        if (a, b, d) == (1, 1, 1 ^ (e * f) ^ (c * g)):
            branches["C"].append(point)
    return branches


def local_triples():
    triples = set()
    for outer in product((-16, 0, 16), repeat=4):
        if sum(outer) != 0:
            continue
        f00, f10, f01, f11 = outer
        numerators = (
            f00 - f10 + f01 - f11,
            f00 + f10 - f01 - f11,
            f00 - f10 - f01 + f11,
        )
        assert all(value % 4 == 0 for value in numerators)
        a, b, c = (value // 4 for value in numerators)
        assert (a + b + c, -a + b - c, a - b - c, -a - b + c) == outer
        triples.add((a, b, c))
    assert len(triples) == 19
    profiles = Counter(tuple(sorted(map(abs, triple))) for triple in triples)
    assert profiles == {(0, 0, 0): 1, (0, 0, 16): 6, (0, 8, 8): 12}
    return sorted(triples)


def walsh_for_mask(character_masks, negative_mask, size):
    return tuple(
        size - 2 * (negative_mask ^ character_mask).bit_count()
        for character_mask in character_masks
    )


def character_masks(points, frequencies):
    return tuple(
        sum(parity(point & frequency) << index for index, point in enumerate(points))
        for frequency in frequencies
    )


def check_graph_saturation(points):
    frequencies = tuple(frequency << 2 for frequency in range(32))
    masks = character_masks(points, frequencies)
    full_mask = (1 << len(points)) - 1
    saturation_masks = {mask ^ offset for mask in masks for offset in (0, full_mask)}
    assert len(saturation_masks) == 64
    saturation_count = Counter()
    for negative_mask in saturation_masks:
        spectrum = walsh_for_mask(masks, negative_mask, len(points))
        saturations = [frequency for frequency, value in zip(frequencies, spectrum)
                       if abs(value) == 16]
        assert len(saturations) == 1
        saturated_layer = (saturations[0] >> 3) & 1
        saturation_count[saturated_layer] += 1
        other_layer = [abs(value)
                       for frequency, value in zip(frequencies, spectrum)
                       if (frequency >> 3) & 1 != saturated_layer]
        assert other_layer == [4] * 16
    assert saturation_count == {0: 32, 1: 32}
    return dict(saturation_count)


def check_two_bit_bent():
    masks = character_masks(tuple(range(4)), tuple(range(4)))
    bent = {mask for mask in range(16) if set(map(abs, walsh_for_mask(masks, mask, 4))) == {2}}
    quadratic = 1 << 3  # ef on (e,f) in little-endian order.
    assert bent == {quadratic ^ mask ^ offset for mask in masks for offset in (0, 15)}
    assert len(bent) == 8
    return len(bent)


def normal_form_phase(point, params):
    a, b, c, d, e, f, g = coordinates(point)
    ac, ag, be, bf, ad, alpha, beta, gamma = params
    if (a, b) == (1, 0):
        return (ac * c) ^ (ad * d) ^ (ag * g) ^ alpha
    if (a, b) == (0, 1):
        return ((1 ^ ad) * d) ^ (be * e) ^ (bf * f) ^ beta
    assert (a, b) == (1, 1)
    return (ad * d) ^ (ac * c) ^ (ag * g) ^ (e * f) ^ (be * e) ^ (bf * f) ^ gamma


def affine_increment(params):
    ac, ag, be, bf, ad, alpha, beta, gamma = params
    a0 = alpha ^ be ^ bf
    b0 = beta ^ ac ^ ag
    # Coordinates are (a,b,c,d,e,f,g,constant).
    coefficients = (b0 ^ gamma, a0 ^ gamma, ac, ad, be, bf, ag, a0 ^ b0 ^ gamma)
    return sum(value << index for index, value in enumerate(coefficients))


def check_normal_form(points, reference):
    masks = character_masks(points, tuple(range(128)))
    full_mask = (1 << len(points)) - 1
    affine_masks = masks + tuple(mask ^ full_mask for mask in masks)
    assert len(set(affine_masks)) == 256
    base = sum(normal_form_phase(point, (0,) * 8) << index
               for index, point in enumerate(points))
    reference_multiplier = affine_masks.index(reference ^ base)
    orbit = {reference ^ affine_mask for affine_mask in affine_masks}
    actual = set()
    increments = set()
    for params in product((0, 1), repeat=8):
        negative_mask = sum(normal_form_phase(point, params) << index
                            for index, point in enumerate(points))
        increment = affine_increment(params)
        assert negative_mask == base ^ affine_masks[increment]
        spectrum = walsh_for_mask(masks, negative_mask, len(points))
        assert set(spectrum) <= {-16, 0, 16} and sum(value != 0 for value in spectrum) == 16
        actual.add(negative_mask)
        increments.add(increment)
    assert increments == set(range(256)) and actual == orbit and len(actual) == 256
    return {
        "parameter_count": len(actual),
        "walsh_valid_count": len(actual),
        "walsh_values_permitted": [-16, 0, 16],
        "walsh_nonzero_per_mask": 16,
        "exact_reference_orbit_match": True,
        "reference_negative_mask": reference,
        "zero_parameter_negative_mask": base,
        "reference_affine_multiplier_little_endian": reference_multiplier,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="R32 exceptional_signed_mask_result.json")
    parser.add_argument("output", type=Path, help="Output audit JSON")
    args = parser.parse_args()
    started = time.perf_counter()
    reference_raw = args.reference.read_bytes()
    reference = json.loads(reference_raw)
    assert len(reference["solutions"]) == 1

    branches = get_branches()
    branch_sizes = {branch: len(points) for branch, points in branches.items()}
    assert branch_sizes == {"A": 8, "B": 8, "C": 16}
    points = sorted(sum(branches.values(), []))
    formula_points = [
        point
        for point in range(128)
        for a, b, c, d, e, f, g in [coordinates(point)]
        if (a * b * d) ^ (a * e * f) ^ (b * c * g)
    ]
    assert points == formula_points == reference["points"]
    triples = local_triples()
    result = {
        "schema": "a28-r034-amplitude-audit-v1",
        "scope": "exact finite checks for the analytic rigidity and phase classification on S in F2^7",
        "reference_sha256": sha256(reference_raw).hexdigest(),
        "coordinates_little_endian": list("abcdefg"),
        "branch_sizes": branch_sizes,
        "outer_frequency_order": [[0, 0], [1, 0], [0, 1], [1, 1]],
        "local_branch_triples": triples,
        "graph_saturation_cases_by_d_frequency": check_graph_saturation(branches["C"]),
        "opposite_d_frequency_absolute_walsh_coefficients": [4] * 16,
        "two_bit_bent_signs": check_two_bit_bent(),
        "normal_form": check_normal_form(points, reference["solutions"][0]),
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "arithmetic": "exact integers, no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    result["elapsed_seconds"] = time.perf_counter() - started
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"verified {len(triples)} local triples, 64 saturated graph characters, and 256 signs")


if __name__ == "__main__":
    main()
