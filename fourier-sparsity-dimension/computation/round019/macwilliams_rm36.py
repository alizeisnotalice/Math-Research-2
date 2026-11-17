#!/usr/bin/env python3
"""Exact independent RM(3,6) weight counts from its RM(2,6) dual."""
import hashlib
import itertools
import json
import math
import platform
import time
from pathlib import Path


def monomial(indices):
    return sum(int(all((x >> j) & 1 for j in indices)) << x for x in range(64))


def rank(vectors):
    p = {}
    for x in vectors:
        while x:
            j = x.bit_length() - 1
            if j in p:
                x ^= p[j]
            else:
                p[j] = x
                break
    return len(p)


def fwht(mask):
    a = [1 - 2 * ((mask >> j) & 1) for j in range(64)]
    h = 1
    while h < 64:
        for i in range(0, 64, 2 * h):
            for j in range(i, i + h):
                l, r = a[j], a[j + h]
                a[j], a[j + h] = l + r, l - r
        h *= 2
    return a


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def main():
    start = time.monotonic()
    b = [monomial(s) for d in range(3) for s in itertools.combinations(range(6), d)]
    c = b + [monomial(s) for s in itertools.combinations(range(6), 3)]
    assert len(b) == rank(b) == 22 and len(c) == rank(c) == 42
    assert all((x & y).bit_count() % 2 == 0 for x in b for y in c)
    quadratics = [monomial(s) for s in itertools.combinations(range(6), 2)]
    dual_hist = [0] * 65
    mask = 0
    for i in range(1 << 15):
        if i:
            mask ^= quadratics[(i & -i).bit_length() - 1]
        for w in fwht(mask):
            dual_hist[(64 - w) // 2] += 1
            dual_hist[(64 + w) // 2] += 1
    assert sum(dual_hist) == 1 << 22
    histogram = []
    for w in range(65):
        numerator = sum(count * sum((-1) ** j * choose(i, j) * choose(64 - i, w - j)
                                    for j in range(w + 1)) for i, count in enumerate(dual_hist))
        assert numerator % (1 << 22) == 0
        count = numerator // (1 << 22)
        assert count >= 0
        histogram.append(count)
    assert sum(histogram) == 1 << 42
    assert histogram == histogram[::-1]
    assert histogram[0] == 1 and not any(histogram[1:8])
    result = {'schema': 'a28-r019-macwilliams-rm36-v1',
              'dual_dimension': 22, 'target_dimension': 42,
              'dual_quadratic_parts_enumerated': 1 << 15,
              'dual_codewords_covered': sum(dual_hist),
              'dual_weight_histogram': dual_hist, 'target_weight_histogram': histogram,
              'target_weight16_count': histogram[16],
              'arithmetic': 'integer FWHT and exact integer Krawtchouk/MacWilliams sums',
              'python': platform.python_version(), 'elapsed_seconds': time.monotonic() - start,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'failed_assertions': 0}
    Path(__file__).with_name('macwilliams-rm36.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['target_weight16_count', 'elapsed_seconds', 'failed_assertions']}))


if __name__ == '__main__':
    main()
