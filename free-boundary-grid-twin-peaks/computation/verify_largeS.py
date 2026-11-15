"""Independent Python-int verification of the n=8 large-S completion bound.

This uses a layer-file enumeration of independent subsets of U, then keeps
only consecutive cardinality layers of the C(I) recurrence in memory.  It is
deliberately independent from the original C++ implementation.
"""
from array import array
from collections import defaultdict
from functools import cache
from math import factorial
from pathlib import Path
import json
import os
import platform
import sys
import time

ROOT = Path(os.environ.get("A022_CERT_DIR", "/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/n8-independent-python"))
ROOT.mkdir(parents=True, exist_ok=True)
LAYER_DIR = ROOT / 'layers'


def graph(n):
    N = n*n
    adj = [0]*N
    for y in range(n):
        for x in range(n):
            v = y*n+x
            for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
                if 0 <= X < n and 0 <= Y < n:
                    adj[v] |= 1 << (Y*n+X)
    closed = [a | (1 << v) for v,a in enumerate(adj)]
    return adj, closed


def connected(mask, adj):
    todo = mask & -mask
    seen = 0
    while todo:
        b = todo & -todo
        todo ^= b
        seen |= b
        v = b.bit_length()-1
        todo |= adj[v] & mask & ~seen
    return seen == mask


def boundary(S, adj):
    b = 0
    q = S
    while q:
        x = q & -q
        q -= x
        b |= adj[x.bit_length()-1]
    return b & ~S


def flush(buffers, handles):
    for k, values in buffers.items():
        if values:
            values.tofile(handles[k])
            values.clear()


def write_independent_layers(U, adj, maxk):
    LAYER_DIR.mkdir(exist_ok=True)
    paths = [LAYER_DIR / f'k{k:02d}.u64' for k in range(maxk+1)]
    for p in paths:
        if p.exists():
            p.unlink()
    handles = [p.open('wb') for p in paths]
    buffers = {k: array('Q') for k in range(maxk+1)}
    count = [0]*(maxk+1)
    vertices = [v for v in range(len(adj)) if U >> v & 1]

    def visit(position, available, chosen, cardinality):
        if position == len(vertices):
            buffers[cardinality].append(chosen)
            count[cardinality] += 1
            if len(buffers[cardinality]) >= 8192:
                buffers[cardinality].tofile(handles[cardinality])
                buffers[cardinality].clear()
            return
        v = vertices[position]
        visit(position+1, available, chosen, cardinality)
        if available >> v & 1:
            visit(position+1, available & ~adj[v] & ~(1 << v), chosen | (1 << v), cardinality+1)

    try:
        visit(0, U, 0, 0)
        flush(buffers, handles)
    finally:
        for f in handles:
            f.close()
    return paths, count


def read_layer(path):
    values = array('Q')
    with path.open('rb') as f:
        values.fromfile(f, os.path.getsize(path)//8)
    return values


def ie_h_and_pair_total(n, S, label):
    adj, closed = graph(n)
    N = n*n
    full = (1 << N)-1
    assert connected(S, adj)
    bnd = boundary(S, adj)
    closedS = S | bnd
    U = full & ~closedS
    R = full & ~S
    assert U & closedS == 0
    assert all((closed[v] & S) == 0 for v in range(N) if U >> v & 1)
    pairs = []
    q = bnd
    while q:
        bw = q & -q
        q -= bw
        w = bw.bit_length()-1
        for y in range(N):
            if adj[w] >> y & 1 and not (closedS >> y & 1):
                pairs.append((w,y,closed[w] | closed[y]))
    assert len(pairs) == 12 if n == 8 else len(pairs) >= 0
    maxk = len([v for v in range(N) if U >> v & 1])
    paths, layer_counts = write_independent_layers(U, adj, maxk)
    fact = factorial(R.bit_count())
    h = fact
    pair_total = fact*len(pairs)
    previous = {0: fact}
    layer_report = [{'k':0, 'independent_sets':1}]
    for k in range(1, maxk+1):
        masks = read_layer(paths[k])
        current = {}
        signed = -1 if k & 1 else 1
        for I in masks:
            neighbourhood = 0
            numerator = 0
            q = I
            while q:
                bit = q & -q
                q -= bit
                v = bit.bit_length()-1
                numerator += previous[I ^ bit]
                neighbourhood |= closed[v]
            assert not (neighbourhood & S)
            C, rem = divmod(numerator, neighbourhood.bit_count())
            assert rem == 0
            current[I] = C
            h += signed*C
            eligible = sum((I & pair_closed) == 0 for _,_,pair_closed in pairs)
            pair_total += signed*C*eligible
        layer_report.append({'k':k, 'independent_sets':len(masks)})
        previous = current
        del masks
    denominator = R.bit_count()*(R.bit_count()-1)
    H, rem = divmod(pair_total, denominator)
    assert rem == 0
    return {
        'label': label, 'n':n, 'S_size':S.bit_count(), 'm':R.bit_count(),
        'U_size':U.bit_count(), 'boundary':bnd.bit_count(), 'dependent_pairs':len(pairs),
        'independent_sets':sum(layer_counts), 'independent_layers':layer_report,
        'hS':h, 'pair_total':pair_total, 'H':H,
    }


def direct_completion_h(n, S):
    adj, _ = graph(n)
    full = (1 << (n*n))-1
    @cache
    def h(T):
        if T == full:
            return 1
        return sum(h(T | (1 << v)) for v in range(n*n)
                   if not (T >> v & 1) and (adj[v] & T))
    return h(S)


def mask(n, cells):
    return sum(1 << (y*n+x) for x,y in cells)


def main():
    # Small independent checks use completion recursion; n=5 leaves only a
    # 4x4 complement, so its memoized state space is at most 2^16.
    small = []
    for n, cells, label in [
        (4, [(0,0),(1,0),(0,1)], 'n4_L3'),
        (5, [(x,0) for x in range(5)] + [(0,y) for y in range(1,5)], 'n5_L9'),
    ]:
        S = mask(n, cells)
        result = ie_h_and_pair_total(n, S, label)
        direct = direct_completion_h(n, S)
        assert result['hS'] == direct
        small.append({'label':label, 'hS':result['hS'], 'direct_hS':direct,
                      'independent_sets':result['independent_sets']})
    # Candidate extended L in the 8x8 grid.
    cells = [(2,0)] + [(x,0) for x in range(3,8)] + [(0,1),(1,1),(2,1)] + [(0,y) for y in range(2,8)]
    S = mask(8, cells)
    result = ie_h_and_pair_total(8, S, 'n8_extended_L')
    gap = result['boundary']*result['H'] - 4*result['hS']
    assert result['S_size'] == 15 and result['boundary'] == 14 and result['m'] == 49
    assert result['U_size'] == 35 and result['independent_sets'] == 3843618
    assert gap > 0
    output = {
        'python':sys.version, 'implementation':'layer-file Python arbitrary integers',
        'small_crosschecks':small, 'n8':result,
        'comparison_gap_boundary_H_minus_4hS':gap,
        'candidate_inequality_boundary_H_le_4hS':False,
    }
    (ROOT / 'largeS-python-result.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'n8_hS':result['hS'], 'n8_H':result['H'], 'gap':gap,
                      'independent_sets':result['independent_sets'], 'verified':True}, indent=2))


if __name__ == '__main__':
    main()
