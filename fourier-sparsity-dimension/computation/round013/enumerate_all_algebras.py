#!/usr/bin/env python3
"""Exact q=4 scan of Boolean subalgebras inside the triangular selector code.

Base C=A=<1>, ell=2.  The code S has functions
h(w0,w1)+z0*a0+z1*a1(w0), hence dimension 7.  For each Boolean
subalgebra E of S containing 1, compute the largest E-submodule contained
in S and test dim(E)+dim(core_E(S)) <= 2Q=8.
"""
from itertools import combinations, product
import argparse,hashlib,json,sys,time
from pathlib import Path
started=time.monotonic()

N=7

def rref_subspaces(n):
    # RREF bases, pivot columns decreasing in integer-bit convention.
    # Generate canonical row spaces by pivot set in ordinary columns 0..n-1.
    for d in range(n+1):
        for piv in combinations(range(n),d):
            free=[(i,j) for i,p in enumerate(piv) for j in range(p+1,n) if j not in piv]
            for vals in range(1<<len(free)):
                rows=[1<<p for p in piv]
                for k,(i,j) in enumerate(free):
                    if vals>>k&1: rows[i]|=1<<j
                yield tuple(rows)

def span(rows):
    out={0}
    for r in rows: out|={x^r for x in tuple(out)}
    return out

def rank(rows):
    b={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in b:x^=b[p]
            else:b[p]=x;break
    return len(b)

# coordinates [1,w0,w1,w0w1,z0,z1,z1w0], domain (z0,z1,w0,w1)
truth=[]
for c in range(N):
    mask=0
    for z0,z1,w0,w1 in product(range(2),repeat=4):
        vals=[1,w0,w1,w0*w1,z0,z1,z1*w0]
        if vals[c]:mask|=1<<(z0|z1<<1|w0<<2|w1<<3)
    truth.append(mask)
def evalcode(v):
    ans=0
    for i in range(N):
        if v>>i&1:ans^=truth[i]
    return ans
lookup={evalcode(v):v for v in range(1<<N)}
def mul(u,v): return lookup[evalcode(u)&evalcode(v)] if (evalcode(u)&evalcode(v)) in lookup else None

best=[];nalgebra=0;nsub=0
for rows in rref_subspaces(N):
    nsub+=1
    E=span(rows)
    if 1 not in E:continue
    if any(mul(x,y) not in E for x in E for y in E):continue
    nalgebra+=1
    # core: all r such that e*r stays in S; linear equations modeled by enumeration
    core=[u for u in range(1<<N) if all(mul(e,u) is not None for e in E)]
    # core set is a vector space; its dimension by finding independent members
    cd=rank(core)
    ed=len(rows)
    best.append((ed+cd,ed,cd,rows))
best.sort(reverse=True)

assert best[0][0]==8
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
result={'status':'passed','all_linear_subspaces':nsub,'boolean_subalgebras':nalgebra,
        'maximum_capacity':best[0][0],'all_maximizers':sum(x[0]==best[0][0] for x in best),
        'maximizer_dimensions':sorted(set((x[1],x[2]) for x in best if x[0]==best[0][0])),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python':sys.version,'elapsed_seconds':time.monotonic()-started,
        'scope':'All unital Boolean subalgebras of the specified seven-dimensional Q4 triangular code; independent from support-partition DP.'}
args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
