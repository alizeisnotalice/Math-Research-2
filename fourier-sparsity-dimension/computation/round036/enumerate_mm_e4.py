"""Enumerate six-variable MM fourth directions over a fixed E4."""
import collections
import hashlib
import itertools
import json
import math
import pathlib
import sys

TABLE = range(64)
ROWS = tuple(range(8))
T_ROWS = (0, 4)
OUTSIDE_ROWS = tuple(x for x in ROWS if x not in T_ROWS)
Q = 8


def truth_mask(predicate):
    return sum(1 << point for point in TABLE if predicate(point & 7, point >> 3))


def hex64(value):
    return f"{value:016x}"


def affine_family(generators):
    result = [0]
    for generator in generators:
        result.extend([value ^ generator for value in result])
    assert len(set(result)) == len(result)
    return sorted(result)


def matrix_rank(vectors):
    basis = {}
    for vector in vectors:
        while vector:
            top = vector.bit_length()
            if top in basis:
                vector ^= basis[top]
            else:
                basis[top] = vector
                break
    return len(basis)


def popcount(word):
    return bin(word).count("1")


def run(output):
    base = truth_mask(lambda x, y: popcount(x & y) % 2)
    block = truth_mask(lambda x, _: x in T_ROWS)
    generators = [
        truth_mask(lambda x, _: x == 0),
        truth_mask(lambda x, _: x == 4),
        truth_mask(lambda x, y: x in T_ROWS and y & 4),
    ]
    e4 = affine_family(generators)
    near_threshold = {}
    permutation_row_changes = collections.Counter()
    all_extension_phases_by_exterior_weight = collections.Counter()

    for pi in itertools.permutations(ROWS):
        # For every E4-coset member to be MM, the y2 slope change on the two
        # indicated rows must swap their two outputs.
        if pi[0] ^ pi[4] != 4:
            continue
        exterior_changed_rows = [x for x in OUTSIDE_ROWS if pi[x] != x]
        exterior_fixed_rows = [x for x in OUTSIDE_ROWS if pi[x] == x]
        moved = len(exterior_changed_rows)
        permutation_row_changes[moved] += 1
        lower = 4 * moved
        for constant_flips in range(len(exterior_fixed_rows) + 1):
            weight = lower + 8 * constant_flips
            all_extension_phases_by_exterior_weight[weight] += (
                (1 << (moved + 2)) * math.comb(len(exterior_fixed_rows), constant_flips)
            )
        if lower > Q:
            continue
        slope = truth_mask(lambda x, y: popcount((pi[x] ^ x) & y) % 2)
        for g in range(1 << len(ROWS)):
            phase = slope ^ truth_mask(lambda x, _: g >> x & 1)
            outside = popcount(phase & ~block)
            if outside > Q:
                continue
            coset = sorted(phase ^ value for value in e4)
            representative = coset[0]
            moving = {pi[0], pi[4]} != set(T_ROWS)
            info = (outside, moving)
            assert near_threshold.get(representative, info) == info
            near_threshold[representative] = info

    rows = []
    for representative, (outside, moving) in sorted(near_threshold.items()):
        if representative == 0:
            continue
        coset = affine_family(generators + [representative])
        assert matrix_rank(generators + [representative]) == 4
        assert outside == Q
        rows.append({
            "representative": hex64(representative),
            "family": [hex64(value) for value in coset],
            "exterior_weight": outside,
            "e4_frequency_support_moves": moving,
            "total_union_weight": popcount(block | representative),
        })
    assert sorted({row["total_union_weight"] for row in rows}) == [24]
    assert len(near_threshold) == 91
    assert all_extension_phases_by_exterior_weight[0] == 8
    assert all_extension_phases_by_exterior_weight[8] == 720
    assert sum(all_extension_phases_by_exterior_weight.values()) == 1474560
    record = {
        "schema": "a28-r036-mm-e4-enumeration-v1",
        "python_version": sys.version,
        "source_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "base_truth_table_hex": hex64(base),
        "e4_physical_support_hex": hex64(block),
        "e4_generators_hex": [hex64(value) for value in generators],
        "e4_family_hex": [hex64(value) for value in e4],
        "phase_count_all_compatible_mm": sum(all_extension_phases_by_exterior_weight.values()),
        "mm_permutation_count": sum(permutation_row_changes.values()),
        "mm_permutations_by_number_of_exterior_row_changes": dict(permutation_row_changes),
        "mm_phases_by_exterior_weight": dict(all_extension_phases_by_exterior_weight),
        "minimum_independent_exterior_weight": min(row["exterior_weight"] for row in rows),
        "near_threshold_independent_coset_count": len(rows),
        "near_threshold_moving_coset_count": sum(row["e4_frequency_support_moves"] for row in rows),
        "near_threshold_stationary_coset_count": sum(not row["e4_frequency_support_moves"] for row in rows),
        "sharp_extensions": rows,
    }
    pathlib.Path(output).write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    print({key: value for key, value in record.items() if key != "sharp_extensions"})


if __name__ == "__main__":
    run(sys.argv[1])
