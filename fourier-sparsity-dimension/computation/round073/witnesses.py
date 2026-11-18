"""Independent exact quotient profiles, simplex symmetry and invariance checks."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations
import json
from pathlib import Path
import time


def character(product):
    return 1 - 2 * (product.bit_count() & 1)


def span(vectors):
    space = {0}
    for v in vectors:
        space |= {x ^ v for x in space}
    return frozenset(space)


@lru_cache(maxsize=None)
def all_subspaces(n):
    ambient = range(1 << n)
    new = {frozenset({0})}
    output = set(new)
    for _ in range(n):
        new = {space | {v ^ x for x in space} for space in new for v in ambient if v not in space}
        output |= new
    return sorted(output, key=lambda S: (len(S), sorted(S)))


def walsh(values, characters, domain):
    return {beta: sum(value * character(beta & x) for value, x in zip(values, domain))
            for beta in characters}


def packed_with_capacity(options, n, capacity):
    @lru_cache(maxsize=None)
    def fit(index, loads):
        if index == len(options):
            return True
        for A in options[index]:
            if all(loads[i] < capacity for i in A):
                after = list(loads)
                for i in A:
                    after[i] += 1
                if fit(index + 1, tuple(after)):
                    return True
        return False

    return fit(0, (0,) * n)


def quotient(table, support, H, record_details=False):
    ambient = range(len(table))
    physical_kernel = [x for x in ambient if all(character(x & a) == 1 for a in H)]
    t_points = sorted({min(x ^ y for y in physical_kernel) for x in ambient})
    labels = sorted({min(a ^ h for h in H) for a in ambient})
    local = [walsh([table[t ^ x] for x in physical_kernel], labels, physical_kernel) for t in t_points]
    L = max(sum(bool(c) for c in row.values()) for row in local)
    M = 0
    options = []
    data = []
    for beta in labels:
        band = sorted(a for a in support if a ^ beta in H)
        if not band:
            assert all(not row[beta] for row in local)
            continue
        W = span(a ^ band[0] for a in band)
        M += len(W)
        classes = {}
        for index, t in enumerate(t_points):
            values = tuple(character(t & w) for w in sorted(W))
            classes.setdefault(values, set()).add(index)
        present = {i for i, row in enumerate(local) if row[beta]}
        available = sorted([tuple(sorted(A)) for A in classes.values() if A & present])
        assert available and all(set(A) <= present for A in available)
        options.append(available)
        if record_details:
            data.append({"label": beta, "band": band, "normal": sorted(W),
                         "available_address_cosets": [[t_points[i] for i in A] for A in available]})
    cap = (sum(len(options[0]) for options in options) + len(t_points) - 1) // len(t_points)
    volume_bound = cap
    options.sort(key=lambda row: (len(row), -len(row[0])))
    while not packed_with_capacity(options, len(t_points), cap):
        cap += 1
    assert cap <= L
    result = {"M": M, "L": L, "packing": cap, "volume_bound": volume_bound,
              "bands": len(options), "raw_product": L * M, "packed_product": cap * M}
    if record_details:
        result["slice_nonzero_labels"] = {t: [beta for beta, coeff in row.items() if coeff]
                                          for t, row in zip(t_points, local)}
        result["band_choices"] = data
    return result


def support_of(table):
    return {point for point, value in walsh(table, range(len(table)), range(len(table))).items() if value}


def orbit_maps():
    maps = []
    for perm in permutations(range(6)):
        mapping = []
        for freq in range(32):
            values = [0, *[(freq >> j) & 1 for j in range(5)]]
            new = sum((values[perm[i]] ^ values[perm[0]]) << (i - 1) for i in range(1, 6))
            mapping.append(new)
        maps.append(mapping)
    return maps


def simplex():
    truth = {0, 1, 2, 4, 8, 16}
    table = [-1 if x in truth else 1 for x in range(32)]
    coefficients = walsh(table, range(32), range(32))
    assert all(w == 32 * (point == 0) - 12 + 4 * point.bit_count()
               for point, w in coefficients.items())
    support = {int(point) for point, coeff in coefficients.items() if coeff}
    assert support == {a for a in range(32) if a.bit_count() != 3}

    catalog = all_subspaces(5)
    profiles = {H: quotient(table, support, H) for H in catalog}
    assert min(row["raw_product"] for row in profiles.values()) == 32
    assert min(row["packed_product"] for row in profiles.values()) == 28

    transforms = orbit_maps()
    todo = set(catalog)
    orbits = []
    for H in catalog:
        if H not in todo:
            continue
        orbit = {frozenset(mapping[x] for x in H) for mapping in transforms}
        assert all(profiles[K] == profiles[H] for K in orbit)
        assert orbit <= todo
        todo -= orbit
        orbits.append({"kernel": sorted(H), "orbit_size": len(orbit), **profiles[H]})
    assert not todo
    assert sum(row['orbit_size'] for row in orbits) == 374
    best_H = next(H for H in catalog if profiles[H]["packed_product"] == 28)

    # Test one ambient GL shear, physical translation and character modulation.
    def transform_frequency(alpha):
        return alpha ^ (16 if alpha & 1 else 0)

    def transpose_input(x):
        return x ^ (1 if x & 16 else 0)

    transformed = [character(13 & x) * table[transpose_input(x) ^ 7] for x in range(32)]
    transformed_support = support_of(transformed)
    assert transformed_support == {transform_frequency(a) ^ 13 for a in support}
    for H in catalog:
        image = frozenset(transform_frequency(a) for a in H)
        assert quotient(transformed, transformed_support, image) == profiles[H]

    # Lift with one unused variable; all kernels, including mixed ones, are tested.
    promoted = table + table
    promoted_support = support_of(promoted)
    assert support == promoted_support
    for H in all_subspaces(6):
        assert quotient(promoted, promoted_support, H) == profiles[frozenset(x for x in H if x < 32)]

    return {"negative_points": sorted(truth), "fourier_coefficients": coefficients,
            "all_split_orbit_representatives": orbits,
            "full_quotient_minima": {"raw": 32, "packed": 28},
            "winning_packed_kernel": sorted(best_H),
            "winning_packed_kernel_details": quotient(table, support, best_H, record_details=True),
            "M26_forced_overlap_details": quotient(table, support, span([3, 5, 9]), record_details=True),
            "affine_transform_kernels_checked": len(catalog),
            "dummy_variable_mixed_kernels_checked": len(all_subspaces(6))}


def examples():
    spectra = []
    for mask in (1, 7, 854):
        table = [1 - 2 * ((mask >> x) & 1) for x in range(16)]
        support = support_of(table)
        assert support == set(range(16))
        spectra.append({"truth_mask": mask, **quotient(table, support, frozenset({0, 1}), True)})
    assert [(r["M"], r["L"], r["packing"]) for r in spectra] == [(16, 8, 7), (16, 8, 4), (16, 4, 4)]

    equal_magnitude_bent = []
    for mask in (854, 857):
        table = [1 - 2 * ((mask >> x) & 1) for x in range(16)]
        assert {abs(c) for c in walsh(table, range(16), range(16)).values()} == {4}
        equal_magnitude_bent.append({"truth_mask": mask,
             **quotient(table, set(range(16)), frozenset({0, 1, 2, 3}), True)})
    assert [(r["L"], r["packing"]) for r in equal_magnitude_bent] == [(1, 1), (4, 1)]

    point_flip = []
    for n in range(2, 6):
        table = [1 - 2 * (x == 0) for x in range(1 << n)]
        support = support_of(table)
        for H in all_subspaces(n):
            a = len(H).bit_length() - 1
            b = n - a
            result = quotient(table, support, H)
            expected = (1 << b) if b >= 2 else 1
            packing = expected - 1 if a and b >= 2 else expected
            assert (result["M"], result["L"], result["packing"]) == (1 << n, expected, packing)
        point_flip.append({"n": n, "all_splits_checked": len(all_subspaces(n))})
    return {"full_support_same_kernel_examples": spectra,
            "equal_magnitude_bent_examples": equal_magnitude_bent,
            "point_flip_all_gl": point_flip}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    results = {"simplex": simplex(), "other_examples": examples()}
    results["elapsed_seconds"] = time.perf_counter() - start
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(results["elapsed_seconds"], len(results["simplex"]["all_split_orbit_representatives"]))
