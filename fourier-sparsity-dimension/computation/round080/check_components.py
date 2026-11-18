#!/usr/bin/env python3
"""Exact Walsh exterior-matroid components in all linear quotients."""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import platform
import random
import struct
import time

SEED = 800828
SAMPLES = {5: 256, 6: 32, 7: 4}
MASK32 = (1 << 32) - 1


def points(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def walsh_mask(n, truth):
    values = [1 - 2 * ((truth >> u) & 1) for u in range(1 << n)]
    stride = 1
    while stride < len(values):
        for start in range(0, len(values), 2 * stride):
            for u in range(start, start + stride):
                a, b = values[u], values[u + stride]
                values[u], values[u + stride] = a + b, a - b
        stride *= 2
    return sum(1 << u for u,v in enumerate(values) if v)


def translate(mask, vector):
    return sum(1 << (u ^ vector) for u in points(mask))


def gaussian(n, k):
    num = den = 1
    for i in range(k):
        num *= (1 << n) - (1 << i)
        den *= (1 << k) - (1 << i)
    return num // den


@lru_cache(maxsize=None)
def spaces(n):
    full = (1 << (1 << n)) - 1
    current = {1}
    for k in range(n + 1):
        if len(current) != gaussian(n,k):
            raise AssertionError(('subspace count',n,k,len(current)))
        yield_values = tuple(sorted(current))
        yield_record = tuple((n-k, H) for H in yield_values)
        if k == 0:
            catalog = list(yield_record)
        else:
            catalog.extend(yield_record)
        following = set()
        for H in current:
            left = full ^ H
            while left:
                u = (left & -left).bit_length() - 1
                coset = translate(H,u)
                following.add(H | coset)
                left ^= coset
        current = following
    return catalog


def reduce_by(value, pivots):
    for i,row in pivots:
        if value & (1 << i):
            value ^= row
    return value


@lru_cache(maxsize=65536)
def h_reducer(H):
    pivots = {}
    for u in points(H):
        v = reduce_by(u,pivots.items())
        if v:
            pivots[v.bit_length()-1] = v
    return tuple(pivots.items())


def compute_components(n, exterior):
    ground = tuple(points(exterior))
    vertex_count = len(ground)
    parent = list(range(vertex_count))
    pivots = [None] * (n+1)

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i,frequency in enumerate(ground):
        column = frequency | (1 << n)
        combination = 0
        while column:
            pivot = column.bit_length() - 1
            existing = pivots[pivot]
            if existing is None:
                pivots[pivot] = (column,combination | (1<<i))
                break
            row, expression = existing
            column ^= row
            combination ^= expression
        if not column:
            for j in points(combination):
                parent[root(j)] = root(i)
    components = {}
    for i,u in enumerate(ground):
        owner = root(i)
        components[owner] = components.get(owner,0) | (1<<u)
    return tuple(sorted(components.values()))


@lru_cache(maxsize=65536)
def small_components(n, exterior):
    return compute_components(n,exterior)


def datasets(n, exhaustive):
    if exhaustive:
        return range(1 << (1 << n))
    rng = random.Random(SEED + n)
    return [rng.getrandbits(1<<n) for _ in range(SAMPLES[n])]


@lru_cache(maxsize=None)
def direction_cut(n):
    return [(1,direction,sum(1<<u for u in range(1<<n) if bin(u & direction).count('1') % 2))
            for direction in range(1,1<<n)]


def audit_support(n,support):
    census = Counter()
    examples = {}
    for t,H in spaces(n):
        exterior = support & ~H
        components = small_components(n,exterior) if n <= 5 else compute_components(n,exterior)
        k = len(components)
        if k > 3 or (t <= 1 and k > 2):
            raise AssertionError((n,t,H,support,components,'too many components'))
        if k == 3:
            reducer = h_reducer(H)
            images = [set(reduce_by(u,reducer) for u in points(c)) for c in components]
            if any(len(span) != 1 for span in images):
                raise AssertionError((n,t,H,support,components,images,'non-singleton quotient component'))
            x,y,z = (next(iter(label)) for label in images)
            if not x or not y or not z or x ^ y ^ z:
                raise AssertionError((n,t,H,support,components,images,'not a 3-circuit of labels'))
            examples.setdefault(t,dict(H_mask=H,component_masks=components,quotient_labels=[x,y,z]))
        census[t,k] += 1
    signature = support & MASK32
    for _,direction,cut in direction_cut(n):
        components = small_components(n,support & cut) if n <= 5 else compute_components(n,support & cut)
        component_count = len(components)
        signature = (signature * 1664525 + component_count + 31*direction) & MASK32
    return census,examples,signature


def audit_dimension(n, exhaustive):
    table_count = 0
    hist = Counter()
    first = {}
    cache = {}
    digest = hashlib.sha256()
    for truth in datasets(n,exhaustive):
        table_count += 1
        support = walsh_mask(n,truth)
        if support not in cache:
            cache[support] = audit_support(n,support)
        census,examples,signature = cache[support]
        hist.update(census)
        for t,example in examples.items():
            first.setdefault(t,dict(truth_mask=truth,support_mask=support,**example))
        digest.update(struct.pack('<I',signature))
    if sum(hist.values()) != table_count*len(spaces(n)):
        raise AssertionError('incorrect number of table-H incidences')
    return dict(n=n,exhaustive=exhaustive,seed=None if exhaustive else SEED+n,tables=table_count,
                subspaces=len(spaces(n)),distinct_supports=len(cache),incidences=table_count*len(spaces(n)),
                n_component_census=[dict(codimension=t,n_components=k,incidences=c)
                                    for (t,k),c in sorted(hist.items())],
                three_components_first_examples=first,
                direction_signature_sha256=digest.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('exhaustive','sample'))
    args = parser.parse_args()
    wall,cpu = time.perf_counter(),time.process_time()
    dimensions = (1,2,3,4) if args.mode == 'exhaustive' else SAMPLES
    records = [audit_dimension(n,args.mode == 'exhaustive') for n in dimensions]
    wall_seconds,cpu_seconds = time.perf_counter()-wall,time.process_time()-cpu
    print(json.dumps(dict(status='pass',mode=args.mode,host=platform.node(),python=platform.python_version(),
                          records=records,wall_seconds=wall_seconds,process_cpu_seconds=cpu_seconds,
                          mean_cpu_percent_of_one_core=100*cpu_seconds/wall_seconds),indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
