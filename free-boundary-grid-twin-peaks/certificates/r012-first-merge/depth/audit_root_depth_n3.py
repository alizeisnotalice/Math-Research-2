#!/usr/bin/env python3
"""Independent 9! permutation audit of the first-merge root-depth histogram."""
import itertools
import json
from collections import Counter
from pathlib import Path

n = 3
N = n*n
FULL = (1 << N) - 1
adj = [0]*N
for y in range(n):
    for x in range(n):
        v=y*n+x
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            X,Y=x+dx,y+dy
            if 0 <= X < n and 0 <= Y < n:
                adj[v] |= 1 << (Y*n+X)

def components(mask):
    out=[]
    while mask:
        comp = mask & -mask
        front = comp
        while front:
            b=front & -front; front -= b
            extra=adj[b.bit_length()-1] & mask & ~comp
            comp |= extra; front |= extra
        out.append(comp)
        mask &= ~comp
    return out

def depth(root, comp):
    x,y=root%n,root//n
    return min(abs(x-v%n)+abs(y-v//n) for v in range(N) if not (comp>>v)&1)

hist=Counter(); two_peak=0
for order in itertools.permutations(range(N)):
    seen=0; births=[]; merged=None
    for v in order:
        before=seen
        old=len(components(before)) if before else 0
        if not (adj[v]&before): births.append(v)
        seen |= 1<<v
        now=len(components(seen))
        if merged is None and old==2 and now==1:
            merged=(seen,v)
    if len(births) != 2:
        continue
    two_peak += 1
    assert merged is not None
    S,w=merged
    A,B=components(S^(1<<w))
    for root in births:
        comp=A if (A>>root)&1 else B
        hist[depth(root,comp)] += 1

got=json.loads(Path('root-depth-n3.json').read_text())['all_roots_histogram']
got={int(k):int(v) for k,v in got.items()}
assert two_peak == 143808
assert dict(hist) == got, (hist,got)
out={'status':'passed','n':3,'two_peak_labelings':two_peak,'root_histogram':dict(sorted(hist.items()))}
Path('audit-root-depth-n3.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
