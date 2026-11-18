"""Exact audit of spectrum losses, root-fold rank, and strict test examples."""

from collections import Counter
from itertools import combinations_with_replacement
import json
import platform


def walsh_from_mask(mask, point_count):
    values = [1 - 2 * ((mask >> point) & 1) for point in range(point_count)]
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for index in range(start, start + width):
                a, b = values[index], values[index + width]
                values[index], values[index + width] = a + b, a - b
        width *= 2
    return values


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            index = row.bit_length() - 1
            if index not in pivots:
                pivots[index] = row
                break
            row ^= pivots[index]
    return len(pivots)


def get_spectrum(mask, point_count):
    values = walsh_from_mask(mask, point_count)
    return {frequency: coefficient for frequency, coefficient in enumerate(values) if coefficient}


def audit_pair(mask0, mask1, point_count, child0=None, child1=None):
    if child0 is None:
        child0 = get_spectrum(mask0, point_count)
        child1 = get_spectrum(mask1, point_count)
    lifted_parent = mask0 | (mask1 << point_count)
    parent = get_spectrum(lifted_parent, 2 * point_count)
    assert {frequency & (point_count - 1) for frequency in parent} == child0.keys() | child1.keys()
    s0, s1, s_parent = len(child0), len(child1), len(parent)
    shared = child0.keys() & child1.keys()
    overlap = len(shared)
    equal_magnitude = sum(abs(child0[beta]) == abs(child1[beta]) for beta in shared)
    assert s_parent == 2 * (s0 + s1) - 2 * overlap - equal_magnitude
    r0 = rank(beta | point_count for beta in child0)
    r1 = rank(beta | point_count for beta in child1)
    projected_rank = rank(beta | point_count for beta in child0.keys() | child1.keys())
    parent_rank = rank(alpha | (2 * point_count) for alpha in parent)
    h = r0 + r1 - projected_rank
    epsilon = parent_rank - projected_rank
    assert h >= 0 and epsilon in (0, 1)
    assert parent_rank + 1 == (r0 + 1) + (r1 + 1) - h - 1 + epsilon
    has_double_fiber = any((beta | point_count) in parent and beta in parent
                           for beta in child0.keys() | child1.keys())
    if has_double_fiber:
        assert epsilon == 1
    return {"s0": s0, "s1": s1, "s_parent": s_parent, "M": overlap,
            "equal_magnitude": equal_magnitude, "h": h, "epsilon": epsilon,
            "rank_discount": h + 1 - epsilon, "has_double_fiber": has_double_fiber}


def mask_from_phase(function, point_count):
    return sum(function(point) << point for point in range(point_count))


def main():
    census = {}
    for dimension in range(4):
        points = 1 << dimension
        number_of_masks = 1 << points
        spectra = [get_spectrum(mask, points) for mask in range(number_of_masks)]
        counts = Counter()
        for mask0, mask1 in combinations_with_replacement(range(number_of_masks), 2):
            result = audit_pair(mask0, mask1, points, spectra[mask0], spectra[mask1])
            counts["all_pairs"] += 1
            counts["strict_gain_from_exact_loss"] += int(result["equal_magnitude"] < result["M"])
            counts["strict_gain_from_fold_rank"] += int(result["epsilon"] == 0)
        census[dimension] = dict(counts)

    # A constant and a two-variable bent child. At C=3 the exact spectral
    # test passes; the coarse 3M spectral loss misses by one squared unit.
    f1 = mask_from_phase(lambda value: int(value == 3), 4)
    gain = audit_pair(0, f1, 4)
    assert gain["s0"] == 1 and gain["s1"] == 4 and gain["h"] == 1
    coarse_s = 2 * (gain["s0"] + gain["s1"]) - 3 * gain["M"]
    assert 9 * coarse_s < (3 * (1 + 2) - 1) ** 2 <= 9 * gain["s_parent"]

    fold = audit_pair(0, 1, 1)
    assert fold["s_parent"] == 1 and fold["h"] == 1 and fold["rank_discount"] == 2
    assert 2 * (1 + 1) - fold["rank_discount"] == 2

    def bit(x, i):
        return (x >> i) & 1

    quadratic = lambda x: (bit(x, 0) * bit(x, 1)) ^ (bit(x, 2) * bit(x, 3))
    left = mask_from_phase(quadratic, 16)
    right = mask_from_phase(lambda x: quadratic(x) ^ (bit(x, 0) * bit(x, 2)), 16)
    nonaffine_graph = audit_pair(left, right, 16)
    assert nonaffine_graph["M"] == nonaffine_graph["equal_magnitude"] == 16
    assert not nonaffine_graph["has_double_fiber"] and nonaffine_graph["epsilon"] == 1

    result = {"python": platform.python_version(), "census_by_quotient_dim": census,
              "strict_exact_loss_example_C3": gain, "strict_root_fold_example_C2": fold,
              "single_fibers_without_root_fold": nonaffine_graph}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
