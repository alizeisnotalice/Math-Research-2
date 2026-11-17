#!/usr/bin/env python3
"""Exact integer mask DP for the row-release fibre identity, n=3..5."""
from __future__ import annotations
from collections import defaultdict
from math import factorial
from pathlib import Path
import json

KNOWN={3:156,4:4213440,5:153450311313600}

def calc(n):
    N=n*n; R=set(range(2*n)); H=tuple(range(2*n,N)); M=len(H); ix={v:i for i,v in enumerate(H)}
    adj=[]
    for v in range(N):
        x,y=v%n,v//n
        adj.append({yy*n+xx for xx,yy in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if 0<=xx<n and 0<=yy<n})
    hnb=[sum(1<<ix[u] for u in adj[v] if u in ix) for v in H]
    rnb=[bool(adj[v]&R) for v in H]
    cur={(0,0):(1,1)} # mask,birth_count -> sigma count, sum product(r_i+i-1)
    for k in range(M):
        nxt=defaultdict(lambda:[0,0])
        for (mask,p),(count,weight) in cur.items():
            for i,v in enumerate(H):
                bit=1<<i
                if mask&bit: continue
                if not (rnb[i] or hnb[i]&mask): continue
                birth=not bool(hnb[i]&mask)
                p2=p+int(birth)
                if birth:
                    assert v//n==2
                    mult=(k+1)+p2-1
                else: mult=1
                z=nxt[(mask|bit,p2)]; z[0]+=count; z[1]+=weight*mult
        cur={key:tuple(value) for key,value in nxt.items()}
    full=(1<<M)-1; byp={p:x for (mask,p),x in cur.items() if mask==full}
    HR=sum(c for c,w in byp.values())
    HT=sum(w*factorial(M+n)//factorial(M+p) for p,(c,w) in byp.items())
    a=(n+1)//2; E=1
    for i in range(1,a+1): E*=2*i-1
    lower=E*factorial(M+n)//factorial(M+a)
    upper=factorial(M+n)//factorial(M+1)
    return {'n':n,'M':M,'H_R':HR,'by_birth_count':{str(p):{'sigma_count':c,'sum_deadline_product':w} for p,(c,w) in byp.items()},
            'H_T_formula':HT,'H_T_reference':KNOWN[n],'formula_matches':HT==KNOWN[n],
            'uniform_lower_per_sigma':lower,'uniform_upper_per_sigma':upper,
            'bounds_hold':HR*lower<=HT<=HR*upper,'max_p_observed':max(byp),'max_p_allowed':a}

out={'scope':'Integer finite verification of an analytic fibre formula; no asymptotic inference.','runs':[calc(n) for n in (3,4,5)]}
Path(__file__).with_name('row-release-dp.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
