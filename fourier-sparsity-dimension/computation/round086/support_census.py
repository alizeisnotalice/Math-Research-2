"""Exact uniform cut capacity on one disjoint shard of the four-bit truth tables."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
import multiprocessing
import platform
import socket
import time


def enlarge(space, x):
    return space.union(y ^ x for y in space)


def subspaces():
    spaces = {frozenset([0])}
    for x in range(1, 32):
        spaces.update(enlarge(flat, x) for flat in tuple(spaces) if x not in flat)
    return [(flat, sum(1 << (x - 16) for x in flat if x >= 16)) for flat in spaces]


@lru_cache(maxsize=None)
def cuts(span):
    return [(mask, (len(span) // len(flat)).bit_length())
            for flat, mask in FLATS if flat < span]


@lru_cache(maxsize=None)
def classify(support):
    span = frozenset([0])
    for i in range(16):
        if support >> i & 1:
            span = enlarge(span, i + 16)
    size = support.bit_count()
    dim_plus_one = len(span).bit_length()
    numerator, denominator = dim_plus_one, size
    for mask, demand in cuts(span):
        outside = (support & ~mask).bit_count()
        if demand * denominator > numerator * outside:
            numerator, denominator = demand, outside
    # Ratio of the minimum norm uniform cover to P; equality characterizes b_+=P.
    return str(Fraction(numerator * size, denominator * dim_plus_one))


def run(interval):
    counts = Counter()
    for table in range(*interval):
        values = [1 - 2 * ((table >> i) & 1) for i in range(16)]
        for half in (1, 2, 4, 8):
            for j in range(16):
                if not j & half:
                    x, y = values[j], values[j + half]
                    values[j], values[j + half] = x + y, x - y
        support = sum(1 << j for j, x in enumerate(values) if x)
        counts[classify(support)] += 1
    return counts


FLATS = subspaces()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('start', 'stop', 'workers'):
        parser.add_argument(name, type=int)
    args = parser.parse_args()
    assert 0 <= args.start < args.stop <= 32768 and args.workers > 0
    started = time.monotonic()
    intervals = [(i, min(i + 128, args.stop)) for i in range(args.start, args.stop, 128)]
    with multiprocessing.get_context('fork').Pool(args.workers) as pool:
        counts = sum(pool.map(run, intervals), Counter())
    assert sum(counts.values()) == args.stop - args.start
    print(json.dumps(dict(start=args.start, stop=args.stop, workers=args.workers,
                          python=platform.python_version(), host=socket.gethostname(),
                          ambient_subspaces=len(FLATS), ratios=counts,
                          seconds=time.monotonic() - started), indent=2))
