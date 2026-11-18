"""Independently check all five-dimensional minimum-word parameter tables."""

from collections import Counter
from itertools import combinations
import json
import platform


def character_mask(label):
    return sum(1 << point for point in range(16) if (label & point).bit_count() % 2)


def main():
    character_masks = [character_mask(label) for label in range(16)]
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a, b in combinations(range(1, 16), 2)})
    planes = [tuple(label for label in range(1, 16) if (label & point).bit_count() % 2 == 0)
              for point in range(1, 16)]
    assert len(lines) == 35 and len(planes) == 15
    result = {}
    for rank in (2, 3, 4):
        number_zero = 1 << (4 - rank)
        zero = (1 << number_zero) - 1
        active = 0xFFFF ^ zero
        counts = Counter()
        checked = 0
        negative = active
        while True:
            coefficients = [
                abs(active.bit_count() - 2 * ((negative ^ character) & active).bit_count())
                for character in character_masks
            ]
            checked += 1
            if (all(value < 8 for value in coefficients[1:])
                and all(sum(coefficients[label] for label in part) < 16 for part in lines)
                and all(sum(coefficients[label] for label in part) < 24 for part in planes)
                and sum(coefficients[1:]) < 40):
                counts[coefficients[0]] += 1
            if negative == 0:
                break
            negative = (negative - 1) & active
        assert checked == 1 << active.bit_count()
        result[rank] = {
            "sign_vectors_checked_including_both_global_signs": checked,
            "constant_absolute_coefficient_census": dict(sorted(counts.items())),
        }
    expected = [{10: 24, 12: 2}, {12: 28, 14: 2}, {13: 30, 15: 2}]
    assert [result[rank]["constant_absolute_coefficient_census"] for rank in (2, 3, 4)] == expected
    print(json.dumps({"python": platform.python_version(), "ranks": result}, indent=2))


if __name__ == "__main__":
    main()
