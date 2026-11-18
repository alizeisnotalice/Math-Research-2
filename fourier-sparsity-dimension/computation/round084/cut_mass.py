"""Enumerate every linear subspace of both actual homogeneous spans."""
from fractions import Fraction
from functools import lru_cache
import json
import multiprocessing as mp
import platform
import random
import sys


@lru_cache(maxsize=None)
def subspaces(n):
    if n == 0:
        return [(0, (0,))]
    previous = subspaces(n - 1)
    answer = list(previous)
    half = 1 << (n - 1)
    for dim, sub in previous:
        unseen = set(range(half))
        while unseen:
            representative = min(unseen)
            coset = tuple(representative ^ v for v in sub)
            unseen.difference_update(coset)
            answer.append((dim + 1, sub + tuple(half + u for u in coset)))
    return answer


def compress(spectrum):
    pivots = {}
    result = {}
    for v, a in spectrum.items():
        word = 0
        while v:
            bit = v.bit_length()
            if bit not in pivots:
                pivots[bit] = (v, word ^ (1 << len(pivots)))
                word = 1 << (len(pivots) - 1)
                break
            pivot, delta = pivots[bit]
            word ^= delta
            v ^= pivot
        result[word] = a
    return result, len(pivots)


def cost(spectrum):
    compressed, n = compress(spectrum)
    assert len(compressed) == len(spectrum)
    weights = [compressed.get(j, 0) for j in range(1 << n)]
    total = sum(weights)
    num, den = 0, 1
    flats = subspaces(n)
    for dim, flat in flats:
        if dim == n:
            continue
        outside = total - sum(weights[v] for v in flat)
        assert outside > 0
        if (n - dim + 1) * den > outside * num:
            num, den = n - dim + 1, outside
    return Fraction(num, den), len(flats) - 1


def trial(seed):
    rng = random.Random(seed)
    tables = subspace_count = strict = 0
    for j in range(48):
        truth = [rng.choice((-1, 1)) for _ in range(8)]
        if j % 4 == 0:
            truth = [1] * 8
        coeffs = {8 | e: abs(sum(y * (-1)**((x & e).bit_count() % 2)
                                  for x, y in enumerate(truth))) for e in range(8)}
        spectrum, n = compress({e: a for e, a in coeffs.items() if a})
        first = True
        while n <= 5:
            k = rng.randrange(1, min(6 - n, n, 2) + 1)
            rows = [rng.randrange(1, 1 << n) for _ in range(k)]
            if k == 2 and rows[0] == rows[1]:
                rows.pop()
            new = {}
            active = False
            for v, mass in spectrum.items():
                b = sum(((v & row).bit_count() % 2) << i for i, row in enumerate(rows))
                if b:
                    active = True
                    for alpha in (0, 1, 2*b, 2*b+1):
                        new[v | (alpha << n)] = mass
                else:
                    new[v] = 2 * mass
            assert active
            old_cost, count = cost(spectrum)
            new_cost, new_count = cost(new)
            # New integer weights use twice the previous common denominator.
            assert 2 * new_cost <= old_cost
            strict += 2 * new_cost < old_cost
            subspace_count += count + new_count
            tables += 1
            spectrum, n = compress(new)
            if first and j % 4 == 0:
                assert 2 * new_cost == old_cost
            first = False
    return dict(seed=seed, lifts=tables, subspaces=subspace_count, strict=strict)


if __name__ == '__main__':
    workers, first = map(int, sys.argv[1:])
    with mp.get_context('fork').Pool(workers) as pool:
        rows = pool.map(trial, range(first, first + workers))
    print(json.dumps(dict(python=platform.python_version(), workers=workers, results=rows), indent=2))
