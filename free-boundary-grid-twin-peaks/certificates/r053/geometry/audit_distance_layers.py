"""Enumerate small distance-layer maps and exact fibers of selected larger outputs."""
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import permutations
from math import comb,perm
from random import Random
import argparse
import json
import os


def grid(n):
    points=range(n*n)
    return tuple(sum(1<<v for v in points if abs(u//n-v//n)+abs(u%n-v%n)==1) for u in points)


def rho(p,q,n):
    return abs(p%n-q%n)+abs(p//n-q//n)


def roots(word,adj):
    seen=0
    peaks=[]
    for v in word:
        if not adj[v]&seen:
            peaks.append(v)
        seen|=1<<v
    return tuple(peaks)


def exact_roots(word,adj,p,q):
    seen=0
    pair=(1<<p)|(1<<q)
    for v in word:
        bit=1<<v
        if bool(adj[v]&seen)==bool(pair&bit):
            return False
        seen|=bit
    return True


def directed(p,q,n):
    left=(p%n,p//n); right=(q%n,q//n)
    reflect=tuple(a>b for a,b in zip(left,right))
    def transform(point):
        return tuple(n-1-a if r else a for a,r in zip(point,reflect))
    return transform(left),transform(right),transform


def move(p,q,n,k,single):
    lo,hi,unreflect=directed(p,q,n)
    caps=(*lo,n-1-hi[0],n-1-hi[1])
    if single:
        assert 2*k<=sum(caps)
        if caps[0]+caps[1]>=k:
            caps=(*caps[:2],0,0)
        else:
            caps=(0,0,*caps[2:])
    parts=[]
    for cap in caps:
        amount=min(cap,k)
        parts.append(amount)
        k-=amount
    assert k==0
    a,b,c,d=parts
    x=(lo[0]-a,lo[1]-b); y=(hi[0]+c,hi[1]+d)
    def path(first,last):
        i,j=first
        route=[first]
        for axis,stop in enumerate(last):
            while (i,j)[axis]!=stop:
                if axis==0:
                    i+=1 if stop>i else -1
                else:
                    j+=1 if stop>j else -1
                route.append((i,j))
        return route
    P,Q=path(x,lo),path(y,hi)
    if not single:
        prefix=[P[0],Q[0],*P[1:],*Q[1:]]
        carry=prefix
    elif len(P)>1:
        prefix=carry=P
    else:
        assert len(Q)>1
        prefix=[P[0],*Q]
        carry=Q
    def index(pair):
        i,j=unreflect(pair)
        return n*j+i
    return tuple(map(index,prefix)),tuple(map(index,carry)),tuple(map(index,(x,y)))


def encode(word,prefix,carry):
    overwritten=set(prefix)
    return prefix+tuple(v for v in word if v not in overwritten)


def cap(n,k,single):
    N=n*n
    if single:
        return (k+1)*(perm(N-1,k)+perm(N-1,k+1))
    return comb(k+3,3)*perm(N-1,k+1)


def iterable_parts(k):
    for a in range(k+1):
        for b in range(k-a+1):
            for c in range(k-a-b+1):
                yield a,b,c,k-a-b-c


def fiber(out,x,y,n,k,single,original):
    adj=grid(n)
    lo,hi,unreflect=directed(x,y,n)
    found=count=candidate_slots=0
    for a,b,c,d in iterable_parts(k):
        i,j=lo[0]+a,lo[1]+b
        I,J=hi[0]-c,hi[1]-d
        if not (0<=i<=I<n and 0<=j<=J<n):
            continue
        ri,rj=unreflect((i,j)); si,sj=unreflect((I,J))
        p,q=n*rj+ri,n*sj+si
        if rho(p,q,n)<2:
            continue
        prefix,P,pair=move(p,q,n,k,single)
        if pair!=(x,y) or out[:len(prefix)]!=prefix:
            continue
        rest=out[len(prefix):]
        others=tuple(v for v in P if v!=p)
        candidate_slots+=perm(n*n-1,len(others))
        for slots in permutations(range(1,n*n),len(others)):
            fixed=dict(zip(slots,others,strict=True))
            tail=iter(rest)
            old=[p]
            for index in range(1,n*n):
                old.append(fixed[index] if index in fixed else next(tail))
            if exact_roots(old,adj,p,q):
                count+=1
                found+=tuple(old)==original
    assert found==1 and count<=candidate_slots<=cap(n,k,single)
    return count,candidate_slots


def exhaustive():
    n=3
    adj=grid(n)
    amounts=Counter()
    histogram=defaultdict(Counter)
    ranked=defaultdict(Counter)
    mapping={}
    for p in range(9):
        for q in range(9):
            d=rho(p,q,n)
            if d<2:
                continue
            mapping[p,q]=[(t,single,move(p,q,n,t-d,single))
                for t in range(d+1,5) for single in (False,True)
                if not single or 2*t<=4+d]
    for word in permutations(range(9)):
        old=roots(word,adj)
        if len(old)!=2:
            continue
        p,q=old
        d=rho(p,q,n)
        amounts[d]+=1
        for t,single,(prefix,P,expected) in mapping[old]:
            out=encode(word,prefix,P)
            assert roots(out,adj)==expected and rho(*expected,n)==t
            histogram[d,t,single][out]+=1
            if not single:
                assert out[:2]==expected
                ranked[d,t,word.index(q)+1][out]+=1
    assert sum(amounts.values())==143808
    rows=[]
    for (d,t,single),groups in sorted(histogram.items()):
        largest=max(groups.values())
        assert largest<=cap(n,t-d,single) and len(groups)<=amounts[t]
        assert sum(groups.values())==amounts[d]
        rows.append(dict(d=d,t=t,single=single,source_count=amounts[d],target_count=amounts[t],
                         distinct_outputs=len(groups),max_true_fiber=largest,proved_cap=cap(n,t-d,single)))
    rank_rows=[]
    for (d,t,rank),h in sorted(ranked.items()):
        upper=comb(t-d+3,3)*perm(7,t-d)
        maximum=max(h.values())
        assert maximum<=upper
        rank_rows.append(dict(d=d,t=t,old_second_peak_rank=rank,max_true_fiber=maximum,proved_cap=upper))
    return rows,rank_rows


def choose_old(n,rng):
    adj=grid(n)
    p=rng.randrange(n*n)
    q=rng.choice([v for v in range(n*n) if rho(p,v,n)>=2])
    word=[p]
    seen=1<<p
    def attach(avoid):
        available=[v for v in range(n*n) if adj[v]&seen and not (1<<v)&(seen|avoid)]
        if not available:
            return None
        return rng.choice(available)
    for _ in range(rng.randrange(12)):
        v=attach((1<<q)|adj[q])
        if v is None:
            break
        word.append(v)
        seen|=1<<v
    word.append(q)
    seen|=1<<q
    while len(word)<n*n:
        v=attach(0)
        assert v is not None
        word.append(v)
        seen|=1<<v
    assert roots(word,adj)==(p,q)
    return tuple(word),p,q


def one_task(arg):
    n,seed,k,single=arg
    rng=Random(seed)
    while True:
        old,p,q=choose_old(n,rng)
        d=rho(p,q,n)
        t=d+k
        D=2*n-2
        if t<=D and (not single or 2*t<=D+d):
            break
    prefix,P,(x,y)=move(p,q,n,k,single)
    out=encode(old,prefix,P)
    assert roots(out,grid(n))==(x,y) and rho(x,y,n)==t
    exact,candidates=fiber(out,x,y,n,k,single,old)
    return dict(n=n,seed=seed,source_rho=d,target_rho=t,single=single,source_low_rank=old.index(q)+1,
                source=list(old),target=list(out),exact_fiber=exact,enumerated_slots=candidates,proved_cap=cap(n,k,single))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=os.cpu_count())
    args=parser.parse_args()
    tasks=[(n,20261118530+100*n+10*k+2*use+j,k,bool(use))
           for n in (4,5,6) for k in (1,2) for use in (0,1) for j in range(2)]
    with ProcessPoolExecutor(max_workers=args.workers) as batch:
        large=list(batch.map(one_task,tasks))
    core,ranked=exhaustive()
    print(json.dumps({'exhaustive_n3':core,'fixed_rank_n3':ranked,'seeded_exact_fibers':large},indent=2))
