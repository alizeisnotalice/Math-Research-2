#!/usr/bin/env python3
"""Exact small regression for affine-mass rectangle halo enclosures.

This is a two-box arithmetic example for the proof interface.  It is not a
complete periodic-halo calculation for a project candidate.
"""
from fractions import Fraction as F
import json
from math import isqrt


def floor_dyadic_root(numerator: F, denominator_power: int) -> tuple[int, F]:
    """Return floor(D*sqrt(numerator)) and its strictly outward D-grid bound."""
    assert numerator > 0
    d = 1 << denominator_power
    # numerator=A/B.  Find l=floor(sqrt(A*d^2/B)) by exact integer squares.
    scaled_num = numerator.numerator * d * d
    scaled_den = numerator.denominator
    lo = isqrt(scaled_num // scaled_den)
    while (lo + 1) ** 2 * scaled_den <= scaled_num:
        lo += 1
    while lo**2 * scaled_den > scaled_num:
        lo -= 1
    assert lo**2 * scaled_den <= scaled_num < (lo + 1) ** 2 * scaled_den
    return lo, F(lo + 1, d)


def union_area(rectangles: list[tuple[F, F, F, F]]) -> F:
    """Exact area of a finite union of axis-parallel rational rectangles."""
    xs = sorted({z for r in rectangles for z in r[:2]})
    ys = sorted({z for r in rectangles for z in r[2:]})
    total = F(0)
    for left, right in zip(xs, xs[1:]):
        for bottom, top in zip(ys, ys[1:]):
            if any(r[0] <= left and right <= r[1] and r[2] <= bottom and top <= r[3]
                   for r in rectangles):
                total += (right - left) * (top - bottom)
    return total


def main() -> None:
    corners = [(F(-1), F(-1)), (F(-1), F(1)), (F(1), F(-1)), (F(1), F(1))]
    masses = [(F(3) + a + b, F(3) - a - b) for a, b in corners]
    assert all(x > 0 and y > 0 and x + y == 6 for x, y in masses)
    max_one = max(x for x, _ in masses)
    max_two = max(y for _, y in masses)
    assert (max_one, max_two) == (F(5), F(5))
    # M1+M2=6 identically, hence their separate maxima cannot be simultaneous.
    assert max_one + max_two != F(6)
    assert not any(x == 5 and y == 5 for x, y in masses)

    # lambda=1, so M/(4 lambda)=5/4 and t=sqrt(5)/2.
    floor_value, q = floor_dyadic_root(F(5, 4), denominator_power=2)
    assert floor_value == 4
    assert q == F(5, 4)
    assert q * q > F(5, 4)  # deliberately strictly outward, including exact roots

    # F=[0,4]^2; point source boxes at (1,1) and (3,1), R=5.
    rectangles = [
        (F(0), F(9, 4), F(0), F(9, 4)),
        (F(7, 4), F(4), F(0), F(9, 4)),
    ]
    enclosed_area = union_area(rectangles)
    assert enclosed_area == 9
    upper_ratio = enclosed_area / 6
    assert upper_ratio == F(3, 2)
    print(json.dumps({
        "status": "passed",
        "corner_masses": [[str(x), str(y)] for x, y in masses],
        "separate_box_maxima": [str(max_one), str(max_two)],
        "no_common_maximizer": True,
        "dyadic_denominator": 4,
        "floor_D_root": floor_value,
        "outward_radius": str(q),
        "union_area": str(enclosed_area),
        "fixed_total_mass": "6",
        "upper_ratio": str(upper_ratio),
        "scope": "two-box envelope arithmetic only; not a complete project halo.",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
