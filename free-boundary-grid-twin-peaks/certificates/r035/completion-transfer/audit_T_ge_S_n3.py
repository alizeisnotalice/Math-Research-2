#!/usr/bin/env python3
"""Exact P3xP3 audit of the transfer lower bound sum T >= sum S."""
from fractions import Fraction
from itertools import permutations
import hashlib, json, platform
from pathlib import Path

HERE=Path(__file__).resolve().parent

def grid(n=3):
    adj=[set() for _ in range(n*n)]
    for y in range(n):
        for x in range(n):
            u=y*n+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                xx,yy=x+dx,y+dy
                if 0<=xx<n and 0<=yy<n:adj[u].add(yy*n+xx)
    return adj

def onepeak(w,adj):
    seen=set()
    for i,v in enumerate(w):
        if i and not (adj[v]&seen):return False
        seen.add(v)
    return True

def stats(w,adj):
    N=len(w); pos=[0]*N
    for rank,v in enumerate(w,1):pos[v]=rank
    root=w[0]
    S=Fraction()
    for v in range(N):
        if v==root or v in adj[root]:continue
        M=N-min(pos[u] for u in adj[v])+1
        S+=sum((Fraction(1,t-1) for t in range(M,N)),Fraction())
    T=0
    for a in range(1,N-1):
        if not (adj[w[a+1]] & set(w[:a])):T+=a
    return S,T

def main():
    adj=grid();ct=0;sumS=Fraction();sumT=0;pointwise_bad=0
    for w in permutations(range(9)):
        if onepeak(w,adj):
            ct+=1;s,t=stats(w,adj);sumS+=s;sumT+=t
            pointwise_bad += s>t
    assert ct==29568 and sumT==79360 and sumT>=sumS
    out={"kind":"exact-P3xP3-total-T-ge-S-transfer-audit","python":platform.python_version(),
         "onepeak_words":ct,"sum_T":sumT,"sum_S":[sumS.numerator,sumS.denominator],
         "sum_T_ge_sum_S":True,"pointwise_S_greater_than_T_words":pointwise_bad,
         "scope":"the target is a total/expectation inequality, not a pointwise claim"}
    text=json.dumps(out,indent=2,sort_keys=True)+'\n';HERE.joinpath('audit-T-ge-S-n3.json').write_text(text)
    print(text,end='');print('output_sha256='+hashlib.sha256(text.encode()).hexdigest())
if __name__=='__main__':main()
