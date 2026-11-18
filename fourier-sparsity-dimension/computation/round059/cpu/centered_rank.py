"""Decide the rank-three centered-indicator obstruction on the four-cube exactly."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import time
from check import affine_flats, binary_rank


def determinant(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    return sum((-1) ** i * value * determinant([row[:i] + row[i + 1:] for row in matrix[1:]])
               for i, value in enumerate(matrix[0]))


@lru_cache(maxsize=None)
def pattern_catalog(seen):
    patterns = tuple(index for index in range(8) if (seen >> index) & 1)
    coordinates = {point: (1, *[(point >> j) & 1 for j in range(3)]) for point in patterns}
    candidate = None
    for basis in combinations(patterns, 4):
        matrix = [coordinates[point] for point in basis]
        det = determinant(matrix)
        if det:
            candidate = basis, matrix, det
            break
    if candidate is None:
        return None
    basis, matrix, det = candidate
    functions = []
    gradients = []
    for values in product(range(2), repeat=4):
        if sum(values) in (0, 4):
            continue
        numerators = [
            determinant([row[:i] + (value,) + row[i + 1:] for row, value in zip(matrix, values)])
            for i in range(4)
        ]
        affine = {point: sum(a * b for a, b in zip(coordinates[point], numerators)) for point in patterns}
        if all(value in (0, det) for value in affine.values()):
            functions.append(tuple(point for point, value in affine.items() if value == det))
            gradients.append(tuple(value - values[0] for value in values[1:]))
    coordinate_masks = [functions.index(tuple(point for point in patterns if (point >> j) & 1))
                        for j in range(3)]
    circuits = []
    for length in range(2, 5):
        for subset in combinations(range(len(functions)), length):
            # All column differences are signed 0/1 vectors. The maximal
            # 3-by-3 binary determinant is two, so primitive minimal
            # nonzero positive dependencies have coefficients in {1,2}.
            for weights in product((1, 2), repeat=length):
                if all(sum(weight * gradients[i][axis] for i, weight in zip(subset, weights)) == 0
                       for axis in range(3)):
                    circuits.append(sum(1 << i for i in subset))
                    break
    target = sum(1 << i for i in coordinate_masks)
    return functions, circuits, target


def audit():
    n = 4
    all_flats = affine_flats(n)
    catalog = {sum(1 << point for point in points): normals for points, normals in all_flats}
    counts = Counter()
    counterexample = None
    full = (1 << (1 << n)) - 1
    proper_flats = [mask for mask in catalog if mask != full]

    for dimension in range(n):
        # All affine flats of a fixed dimension are conjugate by AGL(4,2).
        first = (1 << (1 << dimension)) - 1
        for second, third in combinations((mask for mask in proper_flats if mask != first), 2):
            counts["tested_triples"] += 1
            if binary_rank((*catalog[first], *catalog[second], *catalog[third])) < n:
                continue
            counts["full_normal_rank_triples"] += 1
            # Direct truth lookup: every affine combination is constant
            # on one of these eight sets of addresses.
            atoms = [0] * 8
            for point in range(1 << n):
                mode = sum(((mask >> point) & 1) << j for j, mask in enumerate((first, second, third)))
                atoms[mode] |= 1 << point
            seen = sum(1 << mode for mode, atom in enumerate(atoms) if atom)
            classification = pattern_catalog(seen)
            if classification is None:
                continue
            counts["full_normal_rank_and_real_rank_3_triples"] += 1
            functions, circuits, target = classification
            available = 0
            for index, modes in enumerate(functions):
                mask = sum(atoms[mode] for mode in modes)
                if mask in catalog:
                    available |= 1 << index
            assert target & available == target
            used = 0
            for circuit in circuits:
                if circuit & available == circuit:
                    used |= circuit
            if used & target == target:
                counts["positive_cover_triples"] += 1
                if counterexample is None:
                    counterexample = {"first": first, "second": second, "third": third,
                                      "affine_masks": [sum(atoms[mode] for mode in modes)
                                                       for i, modes in enumerate(functions) if (used >> i) & 1]}
    return {"n": 4, "affine_flat_catalog_size": len(catalog), "counts": counts,
            "counterexample": counterexample,
            "pattern_catalog_cache": pattern_catalog.cache_info()._asdict()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    result = audit()
    result["elapsed_seconds"] = time.perf_counter() - started
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(result)
