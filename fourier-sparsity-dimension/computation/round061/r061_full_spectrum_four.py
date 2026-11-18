"""Integer Walsh verification of the small coefficient cases and the three-point bound."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform


def walsh_masks(n):
    return [sum(1 << x for x in range(1 << n) if (x & xi).bit_count() % 2)
            for xi in range(1 << n)]


def subset_masks(points, size):
    return [sum(1 << x for x in subset) for subset in combinations(points, size)]


def check_five_variables():
    # Affine coordinates put the codimension-two flat P on points 0 through 7.
    masks = walsh_masks(5)
    rows = []
    for r, max_outside in [(5, 0), (7, 2), (9, 4)]:
        for outside_count in range(max(0, r - 8), max_outside + 1):
            inside_masks = subset_masks(range(8), r - outside_count)
            outside_masks = subset_masks(range(8, 32), outside_count)
            norms = Counter()
            for inside in inside_masks:
                for outside in outside_masks:
                    support = inside | outside
                    # The constant frequency must be maximal under the normalization.
                    big = 32 - 2 * r
                    total = big
                    for mask in masks[1:]:
                        z = 2 * abs(r - 2 * (support & mask).bit_count())
                        if z > big:
                            break
                        total += z
                    else:
                        if total < 128:
                            norm = Fraction(total, 32)
                            norms[str(norm)] += 1
                            assert (r, outside_count, norm) in {
                                (5, 0, Fraction(31, 8)), (7, 0, Fraction(29, 8))}
            rows.append({"r": r, "outside": outside_count,
                         "candidates": len(inside_masks) * len(outside_masks),
                         "normalized_below_four": dict(norms)})
    return rows


def enumerated_four_variables():
    counted = Counter()
    masks = walsh_masks(4)
    for support in range(1 << 16):
        r = support.bit_count()
        if r % 2 == 0:
            continue
        z = [abs(16 - 2 * r)]
        z.extend(2 * abs(r - 2 * (support & m).bit_count()) for m in masks[1:])
        if min(z) and sum(z) < 64:
            counted[(str(Fraction(sum(z), 16)), (16 - max(z)) // 2)] += 1
    assert set(counted) == {("11/4", 1), ("13/4", 3), ("7/2", 5)}
    return [{"norm": norm, "distance": distance, "count": count}
            for (norm, distance), count in sorted(counted.items())]


def exact_bounds_at_sixty_four():
    m = 64
    bounds = [2 + Fraction(4*t, t+1) - Fraction(12*t, m) for t in (3, 5, 7)]
    assert min(bounds) > 4
    return list(map(str, bounds))


rows = check_five_variables()
print(json.dumps({
    "python": platform.python_version(),
    "M64_bounds_at_t_3_5_7": exact_bounds_at_sixty_four(),
    "four_variable_odd_supported_masks_with_full_spectrum_below_four": enumerated_four_variables(),
    "five_variable_candidates_after_two_characters_fix_P": rows,
    "five_variable_candidates_checked": sum(row["candidates"] for row in rows),
}, indent=2))
