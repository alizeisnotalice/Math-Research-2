#!/usr/bin/env python3
"""Exact F-partition certificates for a triangle of active E-atoms.

Native bits: x in [0,8), w at bit 3, z at bit 4. No bentness is assumed.
The score is only a certified lower bound on kappa(T); kappa(C) is exact.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import sys
import time

E = (3, 12, 48, 192)
EDGES = ((0, 1), (1, 2), (2, 0), (0, 1), (1, 2))
FULL = (1 << 32) - 1


def basis(vs):
    piv = {}
    for x in vs:
        for j in sorted(piv, reverse=True):
            if x >> j & 1:
                x ^= piv[j]
        if x:
            j = x.bit_length() - 1
            for k in piv:
                if piv[k] >> j & 1:
                    piv[k] ^= x
            piv[j] = x
    return tuple(piv[j] for j in sorted(piv, reverse=True))


def words(b):
    out = [0]
    for g in b:
        out += [x ^ g for x in out]
    return out


def capacity(b):
    code = words(b)
    local_dim = {s: sum(g & ~s == 0 for g in code).bit_length() - 1 for s in code}
    dp = {0: 0}
    pick = {}
    for s in sorted(code, key=int.bit_count):
        if not s:
            continue
        v, a = max((local_dim[a] + 1 + dp[s ^ a], a)
                   for a in code if a & (s & -s) and a & ~s == 0)
        dp[s] = v
        pick[s] = a
    s = 255
    partition = []
    while s:
        a = pick[s]
        partition.append(a)
        s ^= a
    return dp[255], partition


def lift(g, mode):
    return sum(((g >> x) & 1) << (x + 8*w + 16*z)
               for x in range(8) for w in range(2) for z in range(2)
               if mode == 'C' or (mode == 'B' and w) or (mode == 'E' and z))


def candidate(key, cache):
    s = [(key >> x) & 1 for x in range(6)]
    t = [(key >> (6+x)) & 1 for x in range(8)]
    labels = []
    for z in range(2):
        for w in range(2):
            for x in range(8):
                edge = EDGES[x//2] if x < 6 else EDGES[x-3]
                bit = t[x] ^ (z ^ (s[x] & w) if x < 6 else w)
                labels.append(edge[bit])
    cells = [sum(1 << p for p, j in enumerate(labels) if j == v) for v in range(3)]
    assert all(cells) and sum(cells) == FULL
    cc = [sum((labels[x] == j) << x for x in range(8)) for j in range(3)]
    bb = [sum(((labels[x] == j) ^ (labels[x+8] == j)) << x for x in range(8)) for j in range(3)]
    ee = [sum(((labels[x] == j) ^ (labels[x+16] == j)) << x for x in range(8)) for j in range(3)]
    eset = set(words(E))
    assert all(g in eset for g in ee)
    cb = basis([*E, *cc, *(e & b for e in E for b in bb)])
    if cb not in cache:
        cset = set(words(cb))
        bbasis = basis(g for g in cset if all((e & g) in cset for e in E))
        kc, cp = capacity(cb)
        cache[cb] = (bbasis, kc, cp)
    bbase, kc, cp = cache[cb]
    tb = basis([*(lift(g, 'C') for g in cb), *(lift(g, 'B') for g in bbase), *(lift(g, 'E') for g in E)])
    d = len(tb)
    assert all(len(basis([*tb, cell])) == d for cell in cells)
    assert all(cell == lift(c, 'C') ^ lift(b, 'B') ^ lift(e, 'E') for cell, c, b, e in zip(cells, cc, bb, ee))
    # Kernel dimension of restriction to the complement: all supported words,
    # not merely the generators that happen to be supported in a cell.
    local = [d - len(basis(g & ~cell for g in tb)) for cell in cells]
    score = sum(local) + len(cells)
    rec = dict(key=key, C_basis=cb, B_basis=bbase, E_cells=E,
               T_basis=tb, F_cells=cells, F_local_dimensions=local,
               kappa_C=kc, C_optimum_partition=cp, F_score=score,
               gap=score-2*kc, dim_C=len(cb), dim_B=len(bbase), dim_T=d)
    return rec


def block(interval):
    cache = {}
    hist = Counter()
    dims = Counter()
    best = None
    unstable = 0
    for key in range(*interval):
        rec = candidate(key, cache)
        hist[rec['gap']] += 1
        dims[(rec['dim_C'], rec['dim_B'], rec['dim_T'])] += 1
        unstable += rec['dim_B'] < rec['dim_C']
        if best is None or rec['gap'] > best['gap']:
            best = rec
    return dict(checked=interval[1]-interval[0], hist=dict(hist),
                dimensions={str(k): v for k, v in dims.items()}, best=best,
                unstable=unstable)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cases', type=int, default=16384)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    assert 1 <= a.cases <= 16384
    began = time.monotonic()
    hist, dims = Counter(), Counter()
    checked = unstable = 0
    best = None
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(block, (i, min(i+256, a.cases))) for i in range(0, a.cases, 256)]
        for f in as_completed(futures):
            r = f.result()
            hist.update(r['hist']); dims.update(r['dimensions'])
            checked += r['checked']; unstable += r['unstable']
            if best is None or (r['best']['gap'], -r['best']['key']) > (best['gap'], -best['key']):
                best = r['best']
    out = dict(status='complete', checked=checked, workers=a.workers,
               unstable_E_cases=unstable, gap_histogram=dict(hist),
               dimension_histogram=dict(dims), best=best,
               elapsed_seconds=time.monotonic()-began, python=sys.version,
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               scope='Specified 3-cell partitions only, not all algebras of T. Exact base capacities and certified F-module lower bounds. No bentness assertion.')
    a.out.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: out[k] for k in ['status','checked','workers','unstable_E_cases','gap_histogram','elapsed_seconds']}))


if __name__ == '__main__':
    main()
