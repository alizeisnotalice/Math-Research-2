"""Integer checks for maximal bilinear realization on a four-label Boolean circuit."""
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from itertools import product
import json
import multiprocessing as mp
import platform
import sys


def basis(vectors):
    pivots, result = {}, []
    for v in vectors:
        a = v
        while a and a.bit_length() in pivots:
            a ^= pivots[a.bit_length()]
        if a:
            pivots[a.bit_length()] = a
            result.append(v)
    return result


def span(vectors):
    values = {0}
    for v in vectors:
        values |= {x ^ v for x in values}
    return tuple(sorted(values))


@lru_cache(maxsize=None)
def fields(b):
    seen = {(0,)}
    for v in range(1, 1 << b):
        seen.update([span((*i, v)) for i in seen])
    return sorted(seen)


def audit(args):
    b, i, j = args
    options = fields(b)
    checked = attainable = 0
    profiles = Counter()
    for tail in product(options, repeat=2):
        subs = (options[i], options[j], *tail)
        if len(basis(x for line in subs for x in line)) != b:
            continue
        checked += 1
        vectors = [sum(x << (b*j) for j,x in enumerate((*t, t[0]^t[1]^t[2])))
                   for t in product(*subs[:3]) if t[0]^t[1]^t[2] in subs[3]]
        rows = basis(vectors)
        mask = (1 << b) - 1
        projected = [span((v >> (b*j)) & mask for v in rows) for j in range(4)]
        full = list(subs) == projected
        hs = [len(s).bit_length()-1 for s in subs]
        assert len(rows) == sum(hs) - b
        expected = all(len(basis(x for j,s in enumerate(subs) if j != k for x in s)) == b
                       for k in range(4))
        assert full == expected
        if full:
            attainable += 1
            assert b + len(rows) == sum(hs)
            assert Fraction((4 + b + len(rows))**2,sum(4**h for h in hs)) <= 4
            profiles[tuple(sorted(hs))] += 1
    return checked, attainable, profiles


def circuit(k):
    m = 1 << k
    tag = [1 << j for j in range(m-1)] + [(1 << (m-1)) - 1]
    points = []
    for i, t in enumerate(tag):
        if i < k//2:
            addrs = (0, 1 << (2*i), 2 << (2*i))
        else:
            addrs = (m-1 if i == k//2 else 0,)
        points.extend(t*m + a for a in addrs)
    parity = 0
    for p in points:
        parity ^= p
    assert not parity and len(basis(points)) == len(points)-1 == m-1+k
    return len(points)


if __name__ == '__main__':
    jobs, start, stop = map(int,sys.argv[1:])
    tasks = [(3,i,j) for i in range(start,stop) for j in range(16)]
    totals = Counter()
    with mp.get_context('fork').Pool(jobs) as pool:
        records = pool.map(audit,tasks)
    for _,_,p in records:
        totals.update(p)
    print(json.dumps(dict(python=platform.python_version(), workers=jobs, b=3, first=start,
                      stop=stop, spanning_fields=sum(r[0] for r in records),
                      attainable_fields=sum(r[1] for r in records),
                      profiles={','.join(map(str,k)):v for k,v in sorted(totals.items())},
                      circuit_sizes=[circuit(k) for k in (2,4,6,8)]),indent=2))
