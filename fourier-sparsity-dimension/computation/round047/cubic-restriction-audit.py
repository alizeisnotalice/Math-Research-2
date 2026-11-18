"""Check every four-flat in F2^6 and the integer cubic counting bound."""
from collections import Counter
from math import comb, isqrt
import json
import platform


def polynomial_dimension(d):
    return sum(comb(d, t) for t in range(min(d, 3) + 1))


def rank_masks(masks):
    basis = {}
    for mask in masks:
        while mask and mask.bit_length() in basis:
            mask ^= basis[mask.bit_length()]
        if mask:
            basis[mask.bit_length()] = mask
    return len(basis)


def linear_spaces(n, k):
    spaces = {frozenset({0})}
    for _ in range(k):
        extensions = set()
        for space in spaces:
            for v in range(1, 1 << n):
                if v not in space:
                    extensions.add(space.union(x ^ v for x in space))
        spaces = extensions
    return spaces


def affine_four_flats():
    n, k = 6, 4
    monomials = [u for u in range(1 << n) if u.bit_count() <= 3]
    counts = Counter()
    for space in linear_spaces(n, k):
        remaining = set(range(1 << n))
        while remaining:
            origin = min(remaining)
            flat = sorted(origin ^ v for v in space)
            remaining.difference_update(flat)
            masks = (sum(1 << i for i, x in enumerate(flat) if (u & x) == u)
                     for u in monomials)
            counts[rank_masks(masks)] += 1
    assert counts == {polynomial_dimension(k): 2604}
    return counts


def counting_exponent(n):
    # ceil(3*sqrt(n)) using integer arithmetic.
    k = 1 + isqrt(9*n-1)
    assert k <= n
    return (k+1)*n+1-polynomial_dimension(k)


if __name__ == '__main__':
    ranks = affine_four_flats()
    exponents = {n: counting_exponent(n) for n in range(9, 4097)}
    assert max(exponents.values()) < 0
    print(json.dumps({
        "python": platform.python_version(),
        "arithmetic": "exact integer arithmetic",
        "ambient_dimension": 6,
        "flat_dimension": 4,
        "restriction_rank_counts": ranks,
        "counting_exponents_checked_inclusive": [9, 4096],
        "largest_exponent": max(exponents.values()),
        "selected_exponents": {n: exponents[n] for n in (9, 16, 64, 256)}
    }, indent=2))
