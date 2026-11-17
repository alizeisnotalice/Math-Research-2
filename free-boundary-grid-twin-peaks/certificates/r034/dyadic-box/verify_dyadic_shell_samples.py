#!/usr/bin/env python3
"""Exact checks of the shell/reinsertion inequalities on light sampled grid growth words."""
import json, random
from pathlib import Path


def grid(n):
    N=n*n; adj=[set() for _ in range(N)]
    for r in range(n):
        for c in range(n):
            v=r*n+c
            for dr,dc in ((1,0),(-1,0),(0,1),(0,-1)):
                rr,cc=r+dr,c+dc
                if 0<=rr<n and 0<=cc<n: adj[v].add(rr*n+cc)
    return adj


def sample_word(adj, rng):
    r=rng.randrange(len(adj)); seen={r}; w=[r]
    while len(w)<len(adj):
        front=sorted(set().union(*(adj[v] for v in seen))-seen)
        v=front[rng.randrange(len(front))]; seen.add(v); w.append(v)
    return w


def audit_word(n,w,adj):
    N=n*n; root=w[0]; rr,rc=divmod(root,n); pos=[0]*N
    for i,v in enumerate(w):pos[v]=i
    rad=[max(abs(divmod(v,n)[0]-rr),abs(divmod(v,n)[1]-rc)) for v in range(N)]
    D=max(rad); T=[min(pos[v] for v in range(N) if rad[v]==k) for k in range(D+1)]
    records=[]
    for k in range(3,D+1):
        v=w[T[k]]
        early=[u for u in adj[v] if pos[u]<T[k]]
        assert len(early)==1
        p=early[0]
        assert rad[p]==k-1 and p!=root and p not in adj[root]
        lp=min(pos[u] for u in adj[p]); bp=N
        for rank in range(pos[p]+1,N):
            x=w[rank]; earlier=[u for u in adj[x] if pos[u]<rank]
            if earlier==[p]: bp=rank; break
        R=bp-lp-1
        assert lp>=T[k-2] and bp<=T[k] and R<=T[k]-T[k-2]-1
        records.append((k,p,R))
    blocks=[]
    for m in (2,4,8):
        if 2*m<=D:
            block=[R for k,p,R in records if m+1<=k<=2*m]
            assert len(block)==m
            assert sum(block)<=2*T[2*m]<=2*((4*m-1)**2+1)<=32*m*m
            blocks.append(m)
    return D,len(records),blocks


rng=random.Random(20261117); out={"method":"seeded light growth samples; assertions are exact", "samples":{}}
for n,count in ((4,128),(5,256),(6,256),(8,16),(16,8),(32,4)):
    adj=grid(n); Ds=[]; recs=[]; blocks=[]
    for _ in range(count):
        d,c,b=audit_word(n,sample_word(adj,rng),adj); Ds.append(d); recs.append(c); blocks.extend(b)
    out["samples"][str(n)]={"count":count,"D_range":[min(Ds),max(Ds)],
                              "dyadic_blocks_checked":sorted(set(blocks)),
                              "shell_records_total":sum(recs)}
Path(__file__).with_name("dyadic-shell-sample-check.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(json.dumps(out,indent=2,sort_keys=True))
