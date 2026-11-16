#!/usr/bin/env python3
"""Deterministic exact Walsh baseline; standard library only, no float decisions."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
import time


def walsh(v):
    a = list(v)
    h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            for j in range(i, i+h):
                x, y = a[j], a[j+h]
                a[j], a[j+h] = x+y, x-y
        h *= 2
    return a


def basis(vs):
    pivots = {}
    selected = []
    for v in vs:
        w = v
        while w:
            k = w.bit_length()-1
            if k in pivots:
                w ^= pivots[k]
            else:
                pivots[k] = w
                selected.append(v)
                break
    return selected


def verify(v):
    n = len(v)
    a = walsh(v)
    assert walsh(a) == [n*x for x in v]
    assert sum(x*x for x in a) == n*n
    for b in range(1, n):
        assert sum(a[i]*a[i ^ b] for i in range(n)) == 0
    return a


def reduction(v, a):
    bs = basis(i for i,x in enumerate(a) if x)
    h = [None]*(1 << len(bs))
    for x,value in enumerate(v):
        u = sum(((x & b).bit_count() % 2) << j for j,b in enumerate(bs))
        assert h[u] is None or h[u] == value
        h[u] = value
    assert None not in h
    ah = walsh(h)
    assert sum(x != 0 for x in ah) == sum(x != 0 for x in a)
    assert all(ah[1 << j] for j in range(len(bs)))
    f = [(1-x)//2 for x in h]
    af = walsh(f)
    assert all(af[1 << j] for j in range(len(bs)))
    return {'basis':bs, 'reduced_truth':h, 'zero_one_truth':f}


def exhaustive():
    rows = []
    for n in range(5):
        N = 1 << n
        hist = Counter()
        best = {}
        for mask in range(1 << N):
            v = [1-2*((mask >> x) & 1) for x in range(N)]
            a = verify(v)
            s = sum(x != 0 for x in a)
            d = len(basis(i for i,x in enumerate(a) if x))
            hist[d,s] += 1
            # Every truth table is reduced and encoding checked, not just minima.
            cert = reduction(v,a)
            if d not in best or s < best[d]['s']:
                best[d] = {'s':s,'mask':mask,'truth':v,'walsh_numerators':a,**cert}
        rows.append({'n':n,'functions':1 << N,'histogram':[
            {'d':d,'s':s,'count':c} for (d,s),c in sorted(hist.items())],
            'minimum_by_rank':best})
    return rows


def tree_paths(k):
    L = 1 << k
    out = []
    for leaf in range(L):
        node = 0
        path = []
        for level in range(k):
            bit = (leaf >> (k-1-level)) & 1
            path.append((node,bit))
            node = 2*node+1+bit
        out.append(path)
    return out


def tree_spectrum(k):
    L = 1 << k
    spec = {}
    for leaf,path in enumerate(tree_paths(k)):
        for mask in range(1 << k):
            a = 1 << (L-1+leaf)
            sign = 1
            for j,(node,bit) in enumerate(path):
                if (mask >> j) & 1:
                    a ^= 1 << node
                    sign *= 1-2*bit
            assert a not in spec
            spec[a] = sign  # denominator 2**k
    return spec


def tree_checks():
    rows=[]
    for k in range(7):
        L=1 << k
        a=tree_spectrum(k)
        assert len(a)==L*L
        assert len(basis(a))==2*L-1
        conv=Counter()
        for u,x in a.items():
            for v,y in a.items():
                conv[u ^ v]+=x*y
        assert conv[0]==L*L and all(v==0 for u,v in conv.items() if u)
        direct = k <= 3
        if direct:
            paths=tree_paths(k)
            truth=[]
            for x in range(1 << (2*L-1)):
                leaf=next(i for i,path in enumerate(paths)
                          if all(((x >> node)&1)==bit for node,bit in path))
                truth.append(1-2*((x >> (L-1+leaf))&1))
            w=walsh(truth)
            assert {i:v for i,v in enumerate(w) if v} == {
                i:v*(len(truth)//L) for i,v in a.items()}
            reduction(truth,w)
        rows.append({'k':k,'L':L,'s':len(a),'d':2*L-1,
                     'ratio':str(Fraction(len(a),(2*L-1)**2)),
                     'exact_sparse_convolution':True,'truth_table_crosscheck':direct})
    return rows


def support_obstruction():
    rows=[]
    for t in range(2,33):
        support=[(1 << (i+1)) ^ h for i in range(t) for h in (0,1)]
        diffs=Counter(u ^ v for u,v in itertools.combinations(support,2))
        assert min(diffs.values())>=2
        assert len(basis(support))==t+1
        rows.append({'t':t,'s':2*t,'d':t+1,'minimum_difference_pairs':min(diffs.values())})
    return rows


def tree_dp():
    dp=[0,1]
    for L in range(2,257):
        dp.append(min(2*(dp[j]+dp[L-j]) for j in range(1,L)))
    assert all(dp[L]>=L*L for L in range(1,len(dp)))
    assert all(dp[1 << k]==(1 << k)**2 for k in range(9))
    return {'leaf_range':[1,256],'minimum_sparsity':dp[1:]}


def addresses():
    rows=[]
    for k in range(4):
        L=1 << k
        n=k+L
        v=[1-2*((x >> (k+(x & (L-1))))&1) for x in range(1 << n)]
        a=walsh(v)
        assert sum(x!=0 for x in a)==L*L
        assert len(basis(i for i,x in enumerate(a) if x))==k+L
        assert walsh(a)==[(1 << n)*x for x in v]
        reduction(v,a)
        rows.append({'k':k,'s':L*L,'d':k+L,'truth_table_size':len(v)})
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    start=time.monotonic()
    result={'schema':'a28-exact-baseline-v1','python':sys.version,
            'executable':sys.executable,'platform':platform.platform(),
            'arithmetic':'Python arbitrary precision integers; Fraction for ratios',
            'seed':None,'tolerance':0,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'exhaustive':exhaustive(),'trees':tree_checks(),
            'support_only_obstruction':support_obstruction(),'tree_dp':tree_dp()}
    result['addresses']=addresses()
    result['elapsed_seconds']=time.monotonic()-start
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','functions':sum(r['functions'] for r in result['exhaustive']),
                      'minima':[{d:v['s'] for d,v in r['minimum_by_rank'].items()} for r in result['exhaustive']],
                      'trees':result['trees'],'seconds':result['elapsed_seconds']}))


if __name__=='__main__':
    main()
