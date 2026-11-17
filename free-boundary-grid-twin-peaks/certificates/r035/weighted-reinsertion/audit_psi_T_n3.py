#!/usr/bin/env python3
"""Exact P3xP3 cross-check: weighted reinsertion Psi equals T in total."""
from fractions import Fraction
from itertools import permutations
import hashlib, json, platform
from pathlib import Path


HERE=Path(__file__).resolve().parent


def graph(n=3):
    adj=[set() for _ in range(n*n)]
    for y in range(n):
        for x in range(n):
            v=y*n+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                xx,yy=x+dx,y+dy
                if 0<=xx<n and 0<=yy<n:adj[v].add(yy*n+xx)
    return adj


def onepeak(w,adj):
    seen=set()
    for i,v in enumerate(w):
        if i and not (adj[v]&seen):return False
        seen.add(v)
    return True


def psi_T(w,adj):
    N=len(w); pos=[0]*N
    for i,v in enumerate(w,1):pos[v]=i
    root=w[0]; psi=Fraction()
    for p in range(N):
        if p==root or p in adj[root]:continue
        ell=min(pos[u] for u in adj[p]); b=N+1
        for rank in range(pos[p]+1,N+1):
            v=w[rank-1]
            if [u for u in adj[v] if pos[u]<rank]==[p]:
                b=rank;break
        R=b-ell-1
        assert R>=1
        psi+=Fraction(ell-1,R)
    T=0
    for a in range(1,N-1):
        if not (adj[w[a+1]] & set(w[:a])):
            T+=a
    return psi,T


def main():
    adj=graph(); count=0; sumpsi=Fraction();sumT=0; minpsi=None;witness=None
    for w in permutations(range(9)):
        if not onepeak(w,adj):continue
        count+=1; psi,T=psi_T(w,adj);sumpsi+=psi;sumT+=T
        if minpsi is None or psi<minpsi:minpsi,witness=psi,w
    assert count==29568 and sumpsi.denominator==1 and sumpsi==sumT==79360
    out={"kind":"exact-P3xP3-psi-T-total-audit","python":platform.python_version(),
         "onepeak_words":count,"sum_Psi":[sumpsi.numerator,sumpsi.denominator],
         "sum_T":sumT,"min_Psi":[minpsi.numerator,minpsi.denominator],
         "min_witness":list(witness),"Psi_equals_T_total":True,
         "scope":"finite exact cross-check; all-n equality is proved fibrewise in the report"}
    text=json.dumps(out,indent=2,sort_keys=True)+"\n"
    HERE.joinpath("audit-psi-T-n3.json").write_text(text);print(text,end="")
    print("output_sha256="+hashlib.sha256(text.encode()).hexdigest())


if __name__=='__main__':main()
