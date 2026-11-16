#!/usr/bin/env python3
"""Exact integer rank trace for r012 parent-selection equivalence.

It compares the former host CDF target calculation with the new GPU-CDF
algorithm at the integer-operation level.  Mutation follows the Metal
SplitMix64 code verbatim.  This is a small-pool algorithm audit, not an SMC
probability calculation.
"""
import hashlib, json
from pathlib import Path

MASK=(1<<64)-1; C=0x9e3779b97f4a7c15; D=0xd1b54a32d192ed03
K=0xd2b74407b1ce6e93
def mix(x):
    z=(x+C)&MASK; z=((z^(z>>30))*0xbf58476d1ce4e5b9)&MASK
    z=((z^(z>>27))*0x94d049bb133111eb)&MASK; return z^(z>>31)
def nxt(s): return ((s+C)&MASK, mix(s)) # exact old Metal next64: advance then mix(new)
def below(s,b):
    th=(-b)%b
    while True:
        s,x=nxt(s)
        if x>=th:return s,x%b
def key(seed,rank,child,domain): return seed^(rank<<32)^((child*C)&MASK)^((domain*D)&MASK)
def has(m,v): return bool(m[v>>5]&(1<<(v&31)))
def add(m,v):
    m=list(m);m[v>>5]|=1<<(v&31);return tuple(m)
def touch(m,v,n):
    x,y=v%n,v//n
    return (x and has(m,v-1)) or (x+1<n and has(m,v+1)) or (y and has(m,v-n)) or (y+1<n and has(m,v+n))
def allowed(st,n,tau,rank):
    m=st[:5];N=n*n
    return [v for v in range(N) if not has(m,v) and (rank==1 or (not touch(m,v,n) if rank==tau else touch(m,v,n)))]
def choose_state(st,n,tau,rank,child):
    a=allowed(st,n,tau,rank);s=key(SEED,rank,child,2);s,j=below(s,len(a));v=a[j]
    m=add(st[:5],v);r1,r2=st[5:];
    if rank==1:r1=v
    if rank==tau:r2=v
    return (*m,r1,r2)
def run(n,tau,P,seed):
    global SEED;SEED=seed;N=n*n;states=[(0,0,0,0,0,0,0) for _ in range(P)];rows=[]
    for rank in range(1,N+1):
        cnt=[len(allowed(s,n,tau,rank)) for s in states];T=sum(cnt)
        row={'rank':rank,'counts':cnt,'total':T}
        if not T:
            row['dead']=True;rows.append(row);break
        if rank==N:
            moment=0
            for s,w in zip(states,cnt):
                r2=allowed(s,n,tau,rank)[0] if rank==tau else s[6]
                dx=s[5]%n-r2%n;dy=s[5]//n-r2//n;moment+=w*(dx*dx+dy*dy)
            row['terminal_moment']=moment;rows.append(row);break
        cdf=[];q=0
        for w in cnt:q+=w;cdf.append(q)
        old=[];new=[]
        for child in range(P):
            _,u=below(key(seed,rank,child,1),T)
            lo=0;hi=len(cdf)-1
            while lo<hi:
                mid=(lo+hi)//2
                if u<cdf[mid]:hi=mid
                else:lo=mid+1
            old.append(lo)
            # New GPU path: binary-search uploaded 256-particle group CDF,
            # then scan only that group.
            gcdf=[];g=0
            for lo2 in range(0,P,256):
                g+=sum(cnt[lo2:lo2+256]);gcdf.append(g)
            gl=0;gh=len(gcdf)-1
            while gl<gh:
                mid=(gl+gh)//2
                if u<gcdf[mid]:gh=mid
                else:gl=mid+1
            lo2=gl*256;resid=u-(0 if gl==0 else gcdf[gl-1]);acc=0;par=lo2
            for k in range(lo2,min(lo2+256,P)):
                w=cnt[k]
                if resid<acc+w:par=k;break
                acc+=w
            new.append(par)
        assert old==new
        next_old=[choose_state(states[p],n,tau,rank,c) for c,p in enumerate(old)]
        next_new=[choose_state(states[p],n,tau,rank,c) for c,p in enumerate(new)]
        assert next_old==next_new
        row.update(parents=old,states=[list(x) for x in next_old]);rows.append(row);states=next_old
    return rows
cases=[]
for n,taus,P in [(2,range(2,5),513),(3,range(2,10),513),(7,(2,8,47,48),17)]:
    for tau in taus:
        seed=0xA022000000000000+n*1000000+tau*10000
        trace=run(n,tau,P,seed);raw=json.dumps(trace,separators=(',',':')).encode()
        cases.append({'n':n,'tau':tau,'particles':P,'ok':True,'ranks':len(trace),'trace_sha256':hashlib.sha256(raw).hexdigest(),'trace':trace})
out=Path(__file__).with_name('rank-audit-results.json')
out.write_text(json.dumps({'kind':'r012_integer_rank_audit','cases':cases},indent=2)+'\n')
print(out)
