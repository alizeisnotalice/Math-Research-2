"""Enumerate exact cut-mass constants over all four-bit sign tables."""
from collections import Counter
from fractions import Fraction
import json
import multiprocessing as mp
import platform
import sys
import time


def flats():
    # Only affine flats of the sixteen frequency labels can meet the support.
    # Empty is represented by the zero space of homogeneous columns.
    known = {0: 0}
    layer = {1 << i for i in range(16)}
    for rank in range(1, 6):
        following = set()
        for mask in layer:
            known[mask] = rank
            pts = [i for i in range(16) if mask >> i & 1]
            anchor = pts[0]
            for k in range(16):
                if mask >> k & 1:
                    continue
                shifted = sum(1 << (p ^ anchor ^ k) for p in pts)
                following.add(mask | shifted)
        layer = following
    entries = sorted((dim, mask) for mask, dim in known.items())
    return [(d, m & 255, m >> 8) for d, m in entries]


def subset_weights(w):
    vals = [0]
    for x in w:
        vals.extend([a+x for a in vals])
    return vals


def run(args):
    start, end = args
    cuts = flats()
    counts = Counter()
    best = Fraction(0)
    best_table = 0
    begin = time.process_time()
    for table in range(start, end):
        w = [1-2*((table >> j) & 1) for j in range(16)]
        for h in (1, 2, 4, 8):
            for j in range(0, 16, 2*h):
                for i in range(j, j+h):
                    a,b = w[i],w[i+h]
                    w[i],w[i+h] = a+b,a-b
        w = list(map(abs, w))
        lo, hi = subset_weights(w[:8]), subset_weights(w[8:])
        total = sum(w)
        minima = [total] * 6
        for dim, mask_a, mask_b in cuts:
            outside = total - lo[mask_a] - hi[mask_b]
            if outside < minima[dim]:
                minima[dim] = outside
        rank = minima.index(0)
        value = max(Fraction(16*(rank-d+1), minima[d]) for d in range(rank))
        counts[value] += 1
        if value > best:
            best, best_table = value, table
    return counts, best, best_table, time.process_time()-begin


if __name__ == '__main__':
    workers, a, b = map(int,sys.argv[1:])
    tasks = [(j,min(j+128,b)) for j in range(a,b,128)]
    begin = time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        results = pool.map(run,tasks,1)
    elapsed = time.perf_counter()-begin
    hist = sum((x[0] for x in results), Counter())
    best = max(hist)
    print(json.dumps(dict(python=platform.python_version(), workers=workers, start=a,end=b,
                    counts={str(k):v for k,v in sorted(hist.items())},maximum=str(best),
                    first_witness=min(x[2] for x in results if x[1]==best),
                    wall_seconds=elapsed, average_cores=sum(x[3] for x in results)/elapsed),indent=2))
