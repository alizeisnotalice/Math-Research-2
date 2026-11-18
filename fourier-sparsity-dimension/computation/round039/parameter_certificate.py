"""Exact four-parameter tables for the dimension-five proof. Standard library only."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform


def walsh(values):
    out = list(values)
    for width in (1, 2, 4, 8):
        for start in range(0, 16, 2 * width):
            for i in range(start, start + width):
                a, b = out[i], out[i + width]
                out[i], out[i + width] = a + b, a - b
    return out


lines = {tuple(sorted((i, j, i ^ j))) for i in range(1, 16) for j in range(1, i)}
planes = [tuple(i for i in range(1, 16) if (i & w).bit_count() % 2 == 0)
          for w in range(1, 16)]
assert len(lines) == 35 and len(planes) == 15
sidon = (0, 1, 2, 4, 8, 15)
assert {a ^ b for a, b in combinations(sidon, 2)} == set(range(1, 16))


def census(zeros, limits):
    single, line, plane, total = limits
    active = [i for i in range(16) if i not in zeros]
    counts = Counter()
    for mask in range(1 << (len(active) - 1)):
        values = [0] * 16
        for j, x in enumerate(active):
            values[x] = 1 - 2 * ((mask >> j) & 1)
        coeffs = list(map(abs, walsh(values)))
        if max(coeffs[1:]) >= single or sum(coeffs[1:]) >= total:
            continue
        if any(sum(coeffs[i] for i in subset) >= line for subset in lines):
            continue
        if any(sum(coeffs[i] for i in subset) >= plane for subset in planes):
            continue
        counts[coeffs[0]] += 1
    assert counts
    return dict(sorted(counts.items()))


minimum = {}
for n in (1, 2, 4):
    minimum[n] = census(tuple(range(n)), (8, 16, 24, 40))
assert [min(minimum[n]) for n in (1, 2, 4)] == [13, 12, 10]

three_half = {}
for name, zeros in [('point', (0,)), ('line', (0, 1)), ('three', (0, 1, 2)),
                    ('plane', (0, 1, 2, 3)), ('three_times_two', (0, 1, 2, 3, 4, 5)),
                    ('sidon', sidon)]:
    three_half[name] = census(zeros, (8, 14, 20, 32))
assert [min(d) for d in three_half.values()] == [15, 14, 11, 12, 10, 8]

print(json.dumps(dict(python=platform.python_version(), minimum=minimum,
                      three_half=three_half), indent=2))
