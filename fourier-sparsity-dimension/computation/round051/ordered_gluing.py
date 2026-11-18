"""Exact Walsh check of the ordered gluing bound with repeated multi-lift labels."""
from collections import Counter
from fractions import Fraction
import json
from random import Random


def fwht(values):
    values = list(values)
    step = 1
    while step < len(values):
        for i in range(0, len(values), 2 * step):
            for j in range(i, i + step):
                a, b = values[j], values[j + step]
                values[j], values[j + step] = a + b, a - b
        step *= 2
    return values


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            k = x.bit_length()
            if k not in pivots:
                pivots[k] = x
                break
            x ^= pivots[k]
    return len(pivots)


def main():
    rng = Random(5101)
    lines = [(0, 1), (2, 4), (3, 7), (5, 6)]
    output = []
    for k in range(4):
        outer_size = 8 << k
        blocks = outer_size // 2
        data_bits = k + 5
        data_size = 1 << data_bits
        size = outer_size * data_size
        bad_counts = Counter()
        smallest = None
        for sample in range(12):
            truth = [0] * size
            seen = set()
            bad = 0
            for i in range(blocks):
                lo, hi = (x + 8 * (i // 4) for x in lines[i % 4])
                left, right = (0, i + 1) if i < k + 4 else rng.sample(range(data_bits), 2)
                a, b = 1 << left, 1 << right
                signs = [1, 1] if sample == 0 else rng.choices([-1, 1], k=2)
                for addr, freq, sign in zip((lo, hi), (a, b), signs):
                    for y in range(data_size):
                        truth[addr + outer_size * y] = sign * (-1 if freq & y else 1)
                if {a, b} <= seen:
                    bad += 1
                seen.update((a, b))
            assert 0 not in truth and bad == blocks - (k + 4)
            spectrum = fwht(truth)
            D = 1 + rank(size + j for j, x in enumerate(spectrum) if x)
            norm = Fraction(sum(map(abs, spectrum)), size)
            assert D <= 2 * norm + bad - (blocks - 1 - (k + 3))
            bad_counts[bad] += 1
            smallest = 2*norm-D if smallest is None else min(smallest, 2*norm-D)
        output.append(dict(k=k, blocks=blocks, samples=sum(bad_counts.values()),
                           exception_counts=dict(bad_counts), smallest_slack=str(smallest)))
    print(json.dumps(dict(seed=5101, models=output), indent=2))


if __name__ == '__main__':
    main()
