"""Exact checks for the 8-parameter case splits and the six-point projection."""
from collections import Counter
from itertools import combinations, product
import json
import platform


def character(x, y):
    return (-1) ** ((x & y).bit_count() % 2)


def spectrum(values):
    return [sum(a * character(x, y) for x, a in enumerate(values))
            for y in range(len(values))]


def is_line_free(subset):
    return all((a ^ b) not in subset for a, b in combinations(subset, 2))


caps = Counter()
points = set(range(1, 8))
for size in range(8):
    for subset in combinations(points, size):
        cap = set(subset)
        if is_line_free(cap):
            caps[size] += 1
            assert not is_line_free(points - cap)

patterns = {}
for zeros in (1, 2, 3):
    types = Counter()
    for nonzero in product((-1, 1), repeat=8 - zeros):
        for positions in combinations(range(8), zeros):
            signs = iter(nonzero)
            values = [0 if i in positions else next(signs) for i in range(8)]
            walsh = spectrum(values)
            types[tuple(sorted(map(abs, walsh)))] += 1
            if zeros == 2 and all(abs(c) < 4 for c in walsh[1:]):
                assert abs(walsh[0]) == 6
            if zeros == 3:
                assert any(abs(c) >= 3 for c in walsh[1:])
    patterns[zeros] = {str(key): count for key, count in types.items()}

sidon = {0, 4, 8, 13, 14, 15}
assert sorted(a ^ b for a, b in combinations(sidon, 2)) == list(range(1, 16))
kernel = spectrum([int(x in sidon) for x in range(16)])
assert abs(kernel[0]) == 6 and set(map(abs, kernel[1:])) == {2}
print(json.dumps(dict(python=platform.python_version(), cap_counts=caps,
                      walsh_types=patterns, sidon_kernel_l1_numerator=sum(map(abs, kernel)),
                      sidon_kernel_l1_denominator=16), indent=2))
