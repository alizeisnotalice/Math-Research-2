"""Exact four-square lower bound from seam-aligned pairs of local maxima."""
from collections import Counter
from functools import cache
from itertools import permutations
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import sys


def adjacency(n):
    return [
        sum(1 << w for w in range(n*n) if abs(w%n-v%n)+abs(w//n-v//n)==1)
        for v in range(n*n)
    ]


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit
        mask ^= bit


def roots(n):
    N = n*n
    adj = adjacency(n)

    @cache
    def peaks(S):
        if not S:
            return factorial(N)
        closed = S
        total = 0
        for v in bits(S):
            closed |= adj[v.bit_length()-1]
            total += peaks(S ^ v)
        q, remainder = divmod(total, closed.bit_count())
        assert remainder == 0
        return q

    counts = [0] * N
    checked = 0

    def enumerate_sets(remaining, S):
        nonlocal checked
        if S:
            checked += 1
            value = (1 if S.bit_count() % 2 else -1)*peaks(S)
            for v in bits(S):
                counts[v.bit_length()-1] += value
        for v in bits(remaining):
            remaining ^= v
            enumerate_sets(remaining & ~adj[v.bit_length()-1], S | v)

    enumerate_sets((1 << N)-1, 0)
    assert peaks.cache_info().currsize == checked + 1
    return counts, checked + 1


def brute(n):
    N = n*n
    adj = adjacency(n)
    histogram = Counter()
    for word in permutations(range(N)):
        rank = {v:i for i,v in enumerate(word)}
        maxima = [v for v in range(N)
                  if all(rank[v] < rank[w.bit_length()-1] for w in bits(adj[v]))]
        if len(maxima) == 1:
            histogram[maxima[0]] += 1
    return [histogram[v] for v in range(N)]


def main():
    for n in (2,3):
        assert roots(n)[0] == brute(n)
    m = 5
    row, states = roots(m)
    assert all(row[x] == row[m-1-x] == row[(m-1)*m+x] for x in range(m))
    # Allowed seam columns are 0 and 1 on the left and 8 and 9 on the right.
    W = sum(row[x]**2 for x in range(m//2))
    lower = Fraction(W, factorial(m*m)**2)**2
    bound_source = Path(sys.argv[1])
    Z = int(json.loads(bound_source.read_text())['upper_counts'][0]['A2_upper_count'])
    count = lower * factorial(4*m*m)
    assert count.denominator == 1
    result = dict(python=platform.python_version(), n=2*m, min_distance=2*m-2*(m//2)+1,
                  states=states, all_root_counts=row, one_side_intersection_factor=W,
                  constructed_words=count.numerator, unconditional=str(lower),
                  from_upper_sha256=hashlib.sha256(bound_source.read_bytes()).hexdigest(),
                  finite_original_conditional_lower=str(count/Z),
                  decimal=float(count/Z))
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
