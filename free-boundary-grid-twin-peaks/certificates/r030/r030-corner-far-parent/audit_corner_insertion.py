#!/usr/bin/env python3
"""Exhaustive n=3 verification; never an asymptotic claim."""
import hashlib
import itertools as it
import json
import math
import platform
import time
from pathlib import Path


def graph(n):
    return [tuple(w for w in range(n*n)
                  if abs(v//n-w//n)+abs(v % n-w % n)==1)
            for v in range(n*n)]


def peaks(word, adj):
    rank = {v:i for i,v in enumerate(word)}
    return tuple(v for v in word
                 if all(rank[v] < rank[w] for w in adj[v] if w in rank))


def main():
    started=time.monotonic()
    n=3; N=n*n; adj=graph(n); corners=(0,n-1,N-n,N-1)
    family=set(); Bfamily=set(); hole_counts={}; tails={}
    for p in corners:
        a,b=adj[p]
        H=tuple(v for v in range(N) if v not in (p,a,b))
        words=[z for z in it.permutations(H) if len(peaks(z,adj))==1]
        hole_counts[str(p)]=len(words)
        for z in words:
            for slots in it.combinations(range(1,N),3):
                for local in ((p,a,b),(p,b,a)):
                    w=[None]*N; w[0]=z[0]
                    for pos,v in zip(slots,local): w[pos]=v
                    rest=iter(z[1:])
                    for i in range(1,N):
                        if w[i] is None: w[i]=next(rest)
                    w=tuple(w)
                    assert peaks(w,adj)==(z[0],p)
                    assert tuple(v for v in w if v not in (p,a,b))==z
                    assert w not in family
                    family.add(w)
                    if slots[0]==1:
                        Bfamily.add(w)
                        L=min(w.index(a),w.index(b))-1
                        tails[L]=tails.get(L,0)+1
    all_A2=all_B=direct_family=direct_Bfamily=0
    for w in it.permutations(range(N)):
        pk=peaks(w,adj)
        if len(pk)!=2: continue
        all_A2+=1
        is_B=w[1]==pk[1]
        all_B+=is_B
        p=pk[1]
        if p not in corners: continue
        D={p,*adj[p]}
        reduced=tuple(v for v in w if v not in D)
        is_family=(len(peaks(reduced,adj))==1)
        assert is_family==(w in family)
        if is_family:
            direct_family+=1
            direct_Bfamily+=is_B
    h=hole_counts['0']
    assert len(family)==direct_family==8*math.comb(N-1,3)*h
    assert len(Bfamily)==direct_Bfamily==4*(N-2)*(N-3)*h
    assert tails=={t:8*(N-2-t)*h for t in range(1,N-2)}
    assert sum(t*c for t,c in tails.items())==len(family)
    assert all_A2==143808 and all_B==81984
    out={'scope':'exhaustive finite n=3 corner-insertion audit; no asymptotic conclusion',
         'n':n,'one_peak_deleted_corner':hole_counts,'A2':all_A2,'B':all_B,
         'family_A2':len(family),'family_B':len(Bfamily),
         'family_B_L_histogram':tails,'all_images_unique_and_two_peaks':True,
         'full_enumeration_membership_matches':True,
         'python':platform.python_version(),'elapsed_seconds':time.monotonic()-started,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path('corner-insertion-n3.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
