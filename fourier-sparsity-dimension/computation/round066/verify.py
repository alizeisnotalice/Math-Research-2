"""Exact quotient checks for Boolean tables. Run with Python 3.10 or later."""
from collections import Counter
from functools import cache
from itertools import combinations
import json
import platform


def rank(vectors):
    pivots = {}
    for row in vectors:
        while row:
            i = row.bit_length()
            if i not in pivots:
                pivots[i] = row
                break
            row ^= pivots[i]
    return len(pivots)


def support(table, n):
    row = [1 - 2 * ((table >> j) & 1) for j in range(1 << n)]
    width = 1
    while width < len(row):
        for i in range(0, len(row), 2 * width):
            for j in range(i, i + width):
                a, b = row[j], row[j + width]
                row[j], row[j + width] = a + b, a - b
        width *= 2
    return tuple(i for i, x in enumerate(row) if x)


@cache
def relations(labels):
    edges = []
    blocks = [{x} for x in labels]
    for a, b, c in combinations(labels, 3):
        d = a ^ b ^ c
        if d <= c or d not in labels:
            continue
        edge = {a, b, c, d}
        edges.append(sum(1 << x for x in edge))
        connected = [block for block in blocks if block & edge]
        blocks = [block for block in blocks if not block & edge]
        blocks.append(set.union(*connected))
    u, r, count = rank(edges), rank(labels), len(blocks) - 1
    assert r <= count + 2 * u and u <= len(labels) - 1 - r
    return count, u, r


def run(n):
    supports = {support(table, n) for table in range(1 << (1 << n))}
    histogram = Counter()
    for spectrum in supports:
        for mask in range(1 << n):
            bands = Counter(x & mask for x in spectrum)
            q = len(spectrum) - bands[0]
            components, u, r = relations(tuple(sorted(set(bands) | {0})))
            assert components ** 2 <= q and max(0, r - 2 * u) ** 2 <= q
            histogram[u] += 1
    return dict(n=n, distinct_spectra=len(supports), tests=sum(histogram.values()),
                tests_by_relation_rank=dict(sorted(histogram.items())))


if __name__ == '__main__':
    print(json.dumps(dict(python=platform.python_version(),
                          results=[run(n) for n in range(1, 5)]), indent=2))
