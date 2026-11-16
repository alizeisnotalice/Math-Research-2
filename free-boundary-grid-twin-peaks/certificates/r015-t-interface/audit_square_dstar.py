#!/usr/bin/env python3
"""Exact subset-DP audit of the T-shaped full-square D-star obstruction."""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json


def square(n):
    m=(n-1)//2; V={(x,y) for x in range(n) for y in range(n)}
    A={(x,y) for x in range(m) for y in range(n-1)}
    B={(x,y) for x in range(m+1,n) for y in range(n-1)}
    C=tuple(sorted(V-A-B))
    index={v:i for i,v in enumerate(C)}
    adj={v:set() for v in V}
    for x,y in V:
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            q=x+dx,y+dy
            if q in V:adj[x,y].add(q)
    return m,A,B,C,index,adj


def completion(seed,C,index,adj):
    g=len(C); touch=[0]*g
    for i,v in enumerate(C):
        for q in adj[v]:
            if q in index:touch[i]|=1<<index[q]
    initial=sum((1<<i) for i,v in enumerate(C) if adj[v]&seed)
    dp=[0]*(1<<g);dp[0]=1
    for mask in range(1<<g):
        z=dp[mask]
        if not z:continue
        for i in range(g):
            b=1<<i
            if not mask&b and ((touch[i]&mask) or (initial&b)):
                dp[mask|b]+=z
    return dp[-1]


rows=[]
for n in (3,5,7):
    m,A,B,C,index,adj=square(n);g=len(C);a=len(A)
    H_A=completion(A,C,index,adj);H_B=completion(B,C,index,adj)
    W=[v for v in C if adj[v]&A and adj[v]&B]
    F=sum(completion(A|B|{w},tuple(v for v in C if v!=w),
                     {v:i for i,v in enumerate(v for v in C if v!=w)},adj) for w in W)
    # Rebuild helper mappings above rather than rely on C's old bit positions.
    special=(m,n-2)
    Cminus=tuple(v for v in C if v!=special)
    hspecial=completion(A|B|{special},Cminus,{v:i for i,v in enumerate(Cminus)},adj)
    K=comb(2*a+g,a)-comb(a+g-1,a)
    actual=Fraction(comb(2*a,a)*F,H_A*K)
    lower=Fraction(factorial(m+1),g*2**g)
    H_closed=Fraction(2*factorial(g),(m+3)*factorial(m+1))
    F_closed=Fraction((6*m+1)*factorial(g-1),4)
    assert H_A==H_B and H_A<=factorial(g)//factorial(m+1)
    assert hspecial==factorial(g-1) and F>=hspecial
    assert H_A==H_closed and F==F_closed
    assert actual>=lower
    rows.append({'n':n,'m':m,'a':a,'g':g,'H_A':str(H_A),'H_B':str(H_B),
                 'F':str(F),'special_completion':str(hspecial),'K':str(K),
                 'H_closed':str(H_closed),'F_closed':str(F_closed),
                 'F_family_over_Dstar':str(actual),'universal_lower_bound':str(lower)})

# Exact check that the analytic lower bound eventually exceeds one and diverges.
threshold=[]
exact_threshold=[]
for m in range(1,101):
    g=4*m+1; low=Fraction(factorial(m+1),g*2**g)
    if low>1 and not threshold: threshold=[m,str(low)]
    a=2*m*m
    H=Fraction(2*factorial(g),(m+3)*factorial(m+1))
    F=Fraction((6*m+1)*factorial(g-1),4)
    K=comb(2*a+g,a)-comb(a+g-1,a)
    ratio=Fraction(comb(2*a,a))*F/(H*K)
    if ratio>1 and not exact_threshold: exact_threshold=[m,str(ratio)]
out={'status':'passed','direct_subset_DP_rows':rows,
     'lower_bound_first_exceeds_one':threshold,
     'exact_family_over_Dstar_first_exceeds_one':exact_threshold,
     'analytic_check_m_1_to_100':True}
Path('square-dstar-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
