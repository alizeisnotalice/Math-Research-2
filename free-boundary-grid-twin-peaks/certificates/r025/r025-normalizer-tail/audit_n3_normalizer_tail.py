#!/usr/bin/env python3
"""Exact P_3 audit of top-two normalizer versus root-distance bounds."""
from __future__ import annotations
import itertools, json
from fractions import Fraction
from pathlib import Path

n=3; N=n*n
adj=[]
for v in range(N):
    x,y=v%n,v//n; z=[]
    for xx,yy in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
        if 0<=xx<n and 0<=yy<n:z.append(yy*n+xx)
    adj.append(z)

def comps(S):
    todo=set(S); ans=[]
    while todo:
        x=todo.pop(); c={x}; st=[x]
        while st:
            u=st.pop()
            for v in adj[u]:
                if v in todo: todo.remove(v); c.add(v); st.append(v)
        ans.append(c)
    return ans

def births(w):
    seen=set(); ans=[]
    for v in w:
        if not any(x in seen for x in adj[v]): ans.append(v)
        seen.add(v)
    return ans

A=B=0; hist={d:0 for d in range(2,2*n-1)}; invsum=Fraction()
rb_by_d={d:[Fraction(),Fraction()] for d in range(2,2*n-1)}
for w in itertools.permutations(range(N)):
    roots=births(w)
    if len(roots)!=2: continue
    A+=1
    top2=(w[0] in roots and w[1] in roots)
    if top2: B+=1
    seen={w[0]}; first=None
    for v in w[1:]:
        hit=[c for c in comps(seen) if any(x in c for x in adj[v])]
        if len(hit)>=2:
            assert first is None and len(hit)==2
            a,b=map(len,hit); first=(a,b,a+b)
        seen.add(v)
    a,b,q=first
    d=abs(roots[0]%n-roots[1]%n)+abs(roots[0]//n-roots[1]//n)
    m=Fraction(q*(q-1),2*a*b)
    assert d<=q and m>=1 and m>=2*Fraction(d-1,d)
    hist[d]+=1; invsum+=1/m
    if top2:
        # j_i is the one-based first-neighbour rank.  Swapping the labels 1
        # and 2 keeps this pair geometry and swaps the two displayed values.
        j=[]
        for root in w[:2]:
            j.append(next(k+1 for k,v in enumerate(w) if k>=2 and v in adj[root]))
        L=Fraction(j[1]-2)
        Lsym=Fraction(j[0]+j[1]-4,2)
        rb_by_d[d][0]+=L; rb_by_d[d][1]+=Lsym
        if d>2: assert Lsym>=Fraction(3,2)

ratio=Fraction(B,A); einv=invsum/A
dist={d:Fraction(k,A) for d,k in hist.items() if k}
linear=sum((Fraction(d,2*(d-1))*p for d,p in dist.items()),Fraction())
out={'scope':'Exact enumeration P_3; validates formula assumptions on every K=2 word.',
     'A2':A,'B_top2':B,'B_over_A':str(ratio),'mean_inverse_m_under_A2':str(einv),
     'first_merge_identity_holds':einv==ratio,'distance_probability':{str(k):str(v) for k,v in dist.items()},
     'linear_upper_bound':str(linear),'linear_slack':str(linear-ratio),
     'thresholds':[],
     'rootpair_symmetric_Rao_Blackwell_by_distance':
       {str(d):{'sum_L':str(x),'sum_Lsym':str(y),'equal':x==y}
        for d,(x,y) in rb_by_d.items() if x or y}}
for r in range(2,2*n-1):
    tail=sum((p for d,p in dist.items() if d>r),Fraction())
    up=Fraction(2*r,r-1)*(1-ratio)
    point=1-Fraction(r-1,2*r)*tail
    out['thresholds'].append({'r':r,'tail':str(tail),'upper':str(up),'tail_le_upper':tail<=up,
      'E1m_le_pointwise_upper':einv<=point,'mbar':str(Fraction(A,B)),
      'mbar_lower':str(1/(1-Fraction(r-1,2*r)*tail))})
Path(__file__).with_name('n3-normalizer-tail.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
