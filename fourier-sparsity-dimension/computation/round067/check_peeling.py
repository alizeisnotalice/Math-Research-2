"""Exhaustively check the peeling lemma for all four-bit Boolean supports."""
from collections import Counter
from itertools import combinations
import json
import platform


def support(table, size):
    return frozenset(i for i in range(size) if sum(
        1 - 2 * (((table >> j) ^ (i & j).bit_count()) & 1)
        for j in range(size)))


def kernel(labels):
    current = labels
    steps = 0
    while current:
        pairs = list(combinations(current, 2))
        counts = Counter(a ^ b for a, b in pairs)
        new = {x for a, b in pairs if counts[a ^ b] >= 3 for x in (a, b)}
        if new == current:
            break
        current = new
        steps += 1
    return current, steps


if __name__ == '__main__':
    n = 4
    size = 1 << n
    # Build supports by integer Walsh butterflies, using an independent direct
    # transform for one representative per distinct support.
    supports = {}
    for table in range(1 << size):
        row = [1 - 2 * ((table >> j) & 1) for j in range(size)]
        for k in range(n):
            step = 1 << k
            for j in range(size):
                if not j & step:
                    a, b = row[j], row[j | step]
                    row[j], row[j | step] = a+b, a-b
        labels = frozenset(i for i, x in enumerate(row) if x)
        supports.setdefault(labels, table)
    for labels, table in supports.items():
        assert support(table, size) == labels
        retained, _ = kernel(labels)
        if len(labels) > 4:
            assert retained == labels
    # For every possible superset, propagate the union of contained large
    # supports by subset dynamic programming. Check containment after peeling.
    unions = [0] * (1 << size)
    for labels in supports:
        if len(labels) > 4:
            mask = sum(1 << b for b in labels)
            unions[mask] = mask
    hist = Counter()
    for mask in range(1 << size):
        labels = {j for j in range(size) if mask >> j & 1}
        for j in labels:
            unions[mask] |= unions[mask ^ (1 << j)]
        kept, steps = kernel(labels)
        assert unions[mask] & ~sum(1 << j for j in kept) == 0
        hist[steps] += 1
    print(json.dumps(dict(python=platform.python_version(), tables=1 << size,
                          distinct_supports=len(supports), supersets=len(unions),
                          peeling_rounds=dict(sorted(hist.items()))), indent=2))
