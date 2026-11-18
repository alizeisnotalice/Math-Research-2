"""Check minimum and second-weight signed half-band geometry by exact arithmetic."""

from itertools import combinations
import json
from operator import xor
from functools import reduce
import platform


def xor_all(values):
    return reduce(xor, values, 0)


def periods(points, ambient_size):
    point_set = set(points)
    return [direction for direction in range(ambient_size)
            if {direction ^ point for point in points} == point_set]


def exhaustive_four_variables():
    ambient_size, N = 16, 4
    census = {}
    for weight in (4, 5, 6):
        checked = 0
        accepted = 0
        valid_supports = 0
        for support in combinations(range(ambient_size), weight):
            phase_masks = [
                sum(1 << index for index, point in enumerate(support)
                    if (point & frequency).bit_count() % 2)
                for frequency in range(ambient_size)
            ]
            this_support = 0
            for normalized_signs in range(0, 1 << weight, 2):
                checked += 1
                spectrum = []
                for mask in phase_masks:
                    coefficient = weight - 2 * (mask ^ normalized_signs).bit_count()
                    if coefficient not in (0, N, -N):
                        break
                    spectrum.append(coefficient)
                else:
                    accepted += 1
                    this_support += 1
                    dual_support = [frequency for frequency, value in enumerate(spectrum) if value]
                    assert len(dual_support) == weight
                    if weight == N:
                        assert len(periods(support, ambient_size)) == N
                        assert len(periods(dual_support, ambient_size)) == N
                    else:
                        assert weight == 6 and len(periods(support, ambient_size)) == 1
                        assert xor_all(support) == xor_all(dual_support) == 0
                        assert len({a ^ b for a, b in combinations(support, 2)}) == 15
                        assert len({a ^ b for a, b in combinations(dual_support, 2)}) == 15
                        assert normalized_signs.bit_count() % 2
            valid_supports += int(this_support > 0)
        census[weight] = {"sign_assignments_one_global_sign_checked": checked,
                          "accepted_assignments": accepted,
                          "accepted_supports": valid_supports}
    assert [census[k]["accepted_assignments"] for k in (4, 5, 6)] == [560, 0, 7168]
    return census


def sharp_selector_family(r):
    macros, N = 1 << (r - 2), 1 << r
    selector_size = 8 * macros
    exceptional = [(0, 1), (2, 4), (3, 7), (5, 6)]
    aligned = [(point, point ^ 1) for point in range(0, 8, 2)]
    partition = [
        ((address << 3) | a, (address << 3) | b)
        for address in range(macros)
        for a, b in (exceptional if address == 0 else aligned)
    ]
    assert len(partition) == N and sorted(point for pair in partition for point in pair) == list(range(selector_size))
    # Each pair receives a separate data character. A true full spectrum
    # frequency is encoded as (data-label index, selector frequency).
    full_support = {
        (label, frequency)
        for label, (left, right) in enumerate(partition)
        for frequency in range(selector_size)
        if ((left ^ right) & frequency).bit_count() % 2 == 0
    }
    half = {point for point in full_support if point[1] & 1}
    assert len(full_support) == N**2 and len(half) == 3 * N // 2
    different_labels = sum(1 for pair in partition if (pair[0] ^ pair[1]) != 1) * 2
    assert different_labels == 6
    return {"N": N, "full_sparsity": len(full_support), "half_band_sparsity": len(half)}


def main():
    result = {
        "python": platform.python_version(),
        "exhaustive_four_variable_N4": exhaustive_four_variables(),
        "sharp_flat_Boolean_selector_family": [sharp_selector_family(r) for r in range(2, 8)],
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
