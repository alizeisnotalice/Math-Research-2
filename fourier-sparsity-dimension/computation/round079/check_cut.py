#!/usr/bin/env python3
"""Exact affine-flat nullity and ternary-circuit checks."""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import platform
import random
import struct
import time


SEED = 790828
SAMPLE_COUNTS = {5: 512, 6: 32, 7: 4}
MASK32 = (1 << 32) - 1


def parity(x):
    return bin(x).count('1') % 2


def walsh(truth, n):
    values = [1 - 2 * ((truth >> x) & 1) for x in range(1 << n)]
    stride = 1
    while stride < len(values):
        for start in range(0, len(values), 2 * stride):
            for i in range(start, start + stride):
                a, b = values[i], values[i + stride]
                values[i], values[i + stride] = a + b, a - b
        stride *= 2
    return values


def points(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def gf2_basis(columns):
    pivots = {}
    for column in columns:
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]
    return tuple(pivots.values())


def dimension(mask):
    return sum(1 for _ in points(mask)).bit_length() - 1


def translate(mask, direction):
    return sum(1 << (x ^ direction) for x in points(mask))


def gaussian_binomial(n, k):
    numerator = denominator = 1
    for i in range(k):
        numerator *= (1 << n) - (1 << i)
        denominator *= (1 << k) - (1 << i)
    return numerator // denominator


@lru_cache(maxsize=None)
def affine_flats(n):
    size = 1 << n
    full_mask = (1 << size) - 1
    spaces = {1}
    output = []
    for d in range(n + 1):
        if len(spaces) != gaussian_binomial(n, d):
            raise AssertionError((n, d, len(spaces), 'incorrect number of subspaces'))
        flats = set()
        next_spaces = set()
        for space in spaces:
            unvisited = full_mask ^ space
            flats.add(space)
            while unvisited:
                representative = (unvisited & -unvisited).bit_length() - 1
                coset = translate(space, representative)
                flats.add(coset)
                next_spaces.add(space | coset)
                unvisited ^= coset
        codim = n - d
        if len(flats) != gaussian_binomial(n, d) * (1 << codim):
            raise AssertionError((n, d, len(flats), 'incorrect affine-flat count'))
        output.extend((codim, flat) for flat in sorted(flats))
        spaces = next_spaces
    return tuple(output)


@lru_cache(maxsize=262144)
def affine_rank(n, mask):
    return len(gf2_basis((1 << n) | x for x in points(mask)))


def signing_from_table(n, truth):
    return sum(1 << i for i, a in enumerate(walsh(truth, n)) if a)


def dataset(n, exhaustive):
    if exhaustive:
        return range(1 << (1 << n))
    rng = random.Random(SEED + n)
    return [rng.getrandbits(1 << n) for _ in range(SAMPLE_COUNTS[n])]


def evaluate_support(n, support, flats):
    histogram = Counter()
    representatives = {}
    signature = support & MASK32
    for t, flat in flats:
        sample = support & flat
        q = sum(1 for _ in points(sample))
        rank = affine_rank(n, sample)
        nullity = q - rank
        slack = 2 * nullity + (1 << (t + 1)) - q
        if slack < 0:
            raise AssertionError((n, support, t, flat, q, nullity, slack))
        histogram[t, slack] += 1
        if slack == 0:
            representatives.setdefault((t, q, nullity), flat)
        signature = (signature * 1664525 + q + 257 * rank + 65537 * t) & MASK32
    return signature, histogram, representatives


def audit_dimension(n, exhaustive):
    flats = affine_flats(n)
    histogram = Counter()
    equalities = Counter()
    examples = {}
    digest = hashlib.sha256()
    supports = {}
    truth_count = 0
    for truth in dataset(n, exhaustive):
        truth_count += 1
        support = signing_from_table(n, truth)
        if support not in supports:
            supports[support] = evaluate_support(n, support, flats)
        signature, counts, equal = supports[support]
        for key in equal:
            flat = equal[key]
            equalities[key] += 1
            examples.setdefault(key, [truth, flat])
        histogram.update(counts)
        digest.update(struct.pack('<I', signature))
    if sum(histogram.values()) != truth_count * len(flats):
        raise AssertionError('incomplete audit accounting')
    return dict(n=n, exhaustive=exhaustive, tables=truth_count, distinct_supports=len(supports),
                seed=None if exhaustive else SEED+n, affine_flats=len(flats),
                incidences_checked=truth_count*len(flats),
                distinct_support_flat_evaluations=len(supports)*len(flats),
                table_signature_sha256=digest.hexdigest(),
                slack_distribution=[dict(t=t,slack=s,count=count) for (t,s),count in sorted(histogram.items())],
                equality_patterns=[dict(t=t,q=q,nullity=l,tables_with_example=count,
                                        first_truth_and_flat_masks=examples[t,q,l])
                                   for (t,q,l),count in sorted(equalities.items())])


def six_circuit():
    n = 4
    support = (0, 1, 2, 4, 8, 15)
    numerators = (-1, 1, 1, 1, 1, 1)
    if len(gf2_basis((1 << n) | x for x in support)) != len(support) - 1:
        raise AssertionError('incorrect affine nullity')
    for remove in support:
        if len(gf2_basis((1 << n) | x for x in support if x != remove)) != len(support) - 1:
            raise AssertionError('zero-sum set is not a circuit')
    values = [sum(a * (1 - 2 * parity(x & u)) for a,u in zip(numerators,support)) for x in range(1 << n)]
    if set(values) != {-4, 0, 4}:
        raise AssertionError('not ternary')
    g = [a // 4 for a in values]
    f = [a if a else (-1)**z for z in (0,1) for a in g]
    if set(f) != {-1,1}:
        raise AssertionError('not Boolean')
    f_mask = sum((a == -1) << i for i,a in enumerate(f))
    spectrum = walsh(f_mask,n+1)
    selected = spectrum[:(1 << n)]
    expected = dict(zip(support,numerators))
    if any(a != 8 * expected.get(i,0) for i,a in enumerate(selected)):
        raise AssertionError('wrong selected Fourier flat')
    if signing_from_table(n+1,f_mask) & ((1 << (1 << n)) - 1) != sum(1 << x for x in support):
        raise AssertionError('missing circuit support')
    corr = Counter()
    for i,a in enumerate(spectrum):
        if a:
            for j in range(i+1,len(spectrum)):
                if spectrum[j]:
                    corr[i^j] += a*spectrum[j]
    if any(corr.values()) or sum(a*a for a in spectrum) != 32**2:
        raise AssertionError('Boolean convolution failed')
    return dict(n=n+1,t=1,selected_support=support,selected_affine_rank=5,nullity=1,
                q=len(support),bound=2+4,boolean_truth_mask=f_mask,ternary_values=g,
                ternary_value_distribution=dict(Counter(g)),full_integer_walsh=spectrum,
                active_walsh_frequencies=sum(a!=0 for a in spectrum),zero_walsh_convolution_except_origin=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('exhaustive', 'sample', 'circuit'))
    args = parser.parse_args()
    wall, cpu = time.perf_counter(), time.process_time()
    if args.mode == 'circuit':
        data = six_circuit()
    elif args.mode == 'exhaustive':
        data = [audit_dimension(n, True) for n in (1,2,3,4)]
    else:
        data = [audit_dimension(n, False) for n in SAMPLE_COUNTS]
    wall, cpu = time.perf_counter()-wall, time.process_time()-cpu
    print(json.dumps(dict(status='pass',mode=args.mode,records=data,host=platform.node(),python=platform.python_version(),
                          wall_seconds=wall,process_cpu_seconds=cpu,mean_cpu_percent_of_one_core=100*cpu/wall),
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
