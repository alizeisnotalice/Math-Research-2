#!/usr/bin/env python3
"""Exact 3x3 audit of the canonical shortest-path lift inverse DP."""
import itertools
import json
from collections import Counter
from pathlib import Path

n=3; N=n*n
adj=[0]*N
for y in range(n):
    for x in range(n):
        v=y*n+x
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            X,Y=x+dx,y+dy
            if 0<=X<n and 0<=Y<n: adj[v]|=1<<(Y*n+X)
def md(a,b):return abs(a%n-b%n)+abs(a//n-b//n)
def peaks(word):
    pos=[0]*N
    for i,x in enumerate(word):pos[x]=i
    return [x for x in range(N) if not any(pos[y]<pos[x] for y in range(N) if (adj[x]>>y)&1)]
def canonical_path(u,v):
    path=[u]; cur=u
    while cur!=v:
        remain=md(cur,v)
        nxt=min(x for x in range(N) if (adj[cur]>>x)&1 and md(x,v)==remain-1)
        path.append(nxt);cur=nxt
    return tuple(path)
def lift(word):
    ps=peaks(word); assert len(ps)==2
    pos={x:i for i,x in enumerate(word)}
    u,v=sorted(ps,key=pos.get); path=canonical_path(u,v); P=path[1:-1]
    return (u,)+P+tuple(x for x in word if x!=u and x not in set(P)),len(path)-1
def geo_prefix(sigma,d):
    return all((adj[sigma[j]]>>sigma[j+1])&1 for j in range(d-1)) and md(sigma[0],sigma[d-1])==d-1
def inverse_dp(sigma,d):
    assert geo_prefix(sigma,d)
    u=sigma[0]; P=tuple(sigma[1:d]); tau=(u,)+tuple(sigma[d:])
    candidates=[]
    pref=(u,)+P
    for v in range(N):
        if v in pref or md(u,v)!=d or not ((adj[P[-1]]>>v)&1):continue
        if canonical_path(u,v)==pref+(v,):candidates.append(v)
    prefix=[0]
    for x in tau:prefix.append(prefix[-1]|(1<<x))
    inserted=[0]*(1<<len(P))
    for S in range(1,1<<len(P)):
        bit=S&-S; inserted[S]=inserted[S^bit]|(1<<P[bit.bit_length()-1])
    full=(1<<len(P))-1; answer=0
    for v in candidates:
        dp={(1,0):1}
        for _ in range(1,N):
            nd={}
            for (i,S),count in dp.items():
                placed=prefix[i]|inserted[S]
                for j,z in enumerate(P):
                    if S>>j&1 or not (adj[z]&placed):continue
                    key=(i,S|(1<<j));nd[key]=nd.get(key,0)+count
                if i<len(tau):
                    x=tau[i]; born=not (adj[x]&placed)
                    if (x==v and born) or (x!=v and not born):
                        key=(i+1,S);nd[key]=nd.get(key,0)+count
            dp=nd
        answer+=dp.get((len(tau),full),0)
    return answer

actual=Counter(); onepeaks=[]
for w in itertools.permutations(range(N)):
    p=peaks(w)
    if len(p)==1:onepeaks.append(w)
    elif len(p)==2:
        sig,d=lift(w);actual[(sig,d)]+=1
assert len(onepeaks)==29568
rows={};total=Counter()
for sigma in onepeaks:
    for d in range(2,5):
        if not geo_prefix(sigma,d):continue
        exact=inverse_dp(sigma,d); observed=actual[(sigma,d)]
        assert exact==observed,(sigma,d,exact,observed)
        row=rows.setdefault(str(d),{'eligible_outputs':0,'nonzero_outputs':0,'total_fibre':0})
        row['eligible_outputs']+=1;row['total_fibre']+=exact
        if exact:row['nonzero_outputs']+=1
        total[d]+=exact
assert sum(total.values())==143808
out={'n':n,'A1':len(onepeaks),'A2':sum(actual.values()),
 'canonical_rule':'at each step choose the least-numbered neighbor reducing Manhattan distance to v',
 'rows':rows,'total_by_distance':dict(total),'status':'passed'}
Path('exact-inverse-fibre-n3.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
