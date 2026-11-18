"""Audit exact first-contact tables and the rank-preserving path map."""
from collections import Counter,defaultdict,deque
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from hashlib import sha256
from math import perm
from functools import cache
import argparse
import json

from audit_distance_layers import grid,rho,roots

D=Path(__file__).resolve().parent


def tables(paths):
    h=defaultdict(lambda:[0,0])
    for path in paths:
        for row in path.read_text().splitlines():
            j,d,far,count,nonmerge=map(int,row.split())
            assert 2<=j and 2<=d<=far and 0<=nonmerge<=count
            h[j,d,far][0]+=count
            h[j,d,far][1]+=nonmerge
    return dict(h)


def trace(word,q,adj):
    red=tuple(v for v in word if v!=q)
    j0=next(i for i,v in enumerate(red) if adj[v]&(1<<q))
    K=red[:j0]
    n=len(adj)
    U=set(range(n))-set(K)-{v for v in range(n) if any(adj[u]&(1<<v) for u in K)}
    assert q in U
    return red,j0,K,U


def rooted_paths(q,U,adj):
    paths={q:(q,)}
    todo=deque([q])
    while todo:
        v=todo.popleft()
        for u in sorted(U-paths.keys()):
            if adj[v]&(1<<u):
                paths[u]=paths[v]+(u,)
                todo.append(u)
    return paths


def brute(n):
    assert n==3
    adj=grid(n)
    result=defaultdict(lambda:[0,0])
    stats=Counter()
    balance=Counter()
    for word in permutations(range(n*n)):
        rs=roots(word,adj)
        if len(rs)!=2:
            continue
        p,q=rs
        red,j0,K,U=trace(word,q,adj)
        next_vertex=red[j0]
        paths=rooted_paths(q,U,adj)
        far=max(paths,key=lambda v:(rho(p,v,n),v))
        d=rho(p,q,n)
        key=j0+1,d,rho(p,far,n)
        result[key][0]+=1
        if next_vertex in U:
            swap=tuple(next_vertex if v==q else q if v==next_vertex else v for v in word)
            assert roots(swap,adj)==(p,next_vertex)
            restored=trace(swap,next_vertex,adj)
            assert restored[1:3]==(j0,K) and restored[0][j0]==q
            result[key][1]+=1
            balance[K,q,next_vertex]+=1
        if far!=q:
            route=paths[far]
            s=set(route)
            new_red=K+tuple(reversed(route[:-1]))+tuple(v for v in red[j0:] if v not in s)
            pos=word.index(q)
            out=new_red[:pos]+(far,)+new_red[pos:]
            assert roots(out,adj)==(p,far)
            assert trace(out,far,adj)[1:3]==(j0,K)
            slots={i:v for i,v in enumerate(red[j0:]) if v in s}
            rest=iter(new_red[j0+len(route)-1:])
            rebuilt=K+tuple(slots[i] if i in slots else next(rest) for i in range(len(red)-j0))
            assert rebuilt==red
            stats['nontrivial_path_moves']+=1
    assert all(balance[K,y,q]==value for (K,q,y),value in balance.items())
    stats['oriented_reversible_edge_counts']=len(balance)
    return dict(result),stats



def all_path_bound():
    adj=grid(3)
    full=(1<<len(adj))-1
    nodes=lambda s: [v for v in range(len(adj)) if s&(1<<v)]
    def ring(s):
        union=0
        for v in nodes(s):
            union|=adj[v]
        return union
    @cache
    def suffix(s):
        return 1 if s==full else sum(suffix(s|(1<<v)) for v in nodes(ring(s)&~s))
    cases=0
    for K in range(1,full+1):
        reached=K&-K
        while extra:=ring(reached)&K&~reached:
            reached|=extra
        if reached!=K:
            continue
        U=full&~(K|ring(K))
        if not U:
            continue
        upper=9-K.bit_count()-1
        H={q:sum(suffix(K|(1<<q)|(1<<u)) for u in nodes(adj[q])) for q in nodes(U)}
        for y in H:
            def layer(q,S,k):
                if not k:
                    return H[q]
                return sum(layer(u,S|(1<<u),k-1) for u in nodes(adj[q]&U&~S))
            for k in range(1, U.bit_count()):
                sum_all=layer(y,1<<y,k)
                assert sum_all<=perm(upper,k)*H[y]
                cases+=sum_all>0
    return cases

def summary(n,actual,exact):
    counts=Counter()
    distant=[]
    total=sum(r[0] for r in actual.values())
    assert total==int(exact['A2'])
    for (j,d,f),(count,v) in actual.items():
        counts[d]+=count
    assert counts=={int(d):int(w) for d,w in exact['distance_counts'].items() if int(w)}
    for target in sorted({max(3,n//2+1),n-1,2*n-2}):
        before=sum(c for (j,d,f),(c,_) in actual.items() if d<target)
        if before==0:
            continue
        reach=sum(c for (j,d,f),(c,_) in actual.items() if d<target<=f)
        distant.append({'target_at_least':target,'source_below_target':before,'reachable':reach,
                        'reachable_fraction':str(Fraction(reach,before))})
    by_rank=[]
    for j in sorted({k[0] for k in actual}):
        for target in (max(3,n//2+1),n-1):
            m=sum(c for (i,d,f),(c,_) in actual.items() if i==j and d<target)
            if m:
                r=sum(c for (i,d,f),(c,_) in actual.items() if i==j and d<target<=f)
                by_rank.append({'rank':j,'target_at_least':target,'source':m,'reachable':r})
    return {'n':n,'cells':len(actual),'A2':total,'nonmerge_first_contact':sum(r[1] for r in actual.values()),
            'reachability':distant,'by_rank':by_rank}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('old_exact',type=Path)
    p.add_argument('old_n5',type=Path)
    args=p.parse_args()
    old={r['n']:r for r in json.loads(args.old_exact.read_text())['results']}
    old[5]=json.loads(args.old_n5.read_text())
    actual={n:tables([D/f'contact-n{n}.tsv'] if n<5 else [D/'n5-local.tsv',D/'n5-edge.tsv']) for n in (3,4,5)}
    for name, data in (('n5-edge','n5-edge.tsv'),('n5-local','n5-local.tsv'),('n5-mac','n5-edge.tsv')):
        receipt=json.loads((D/f'{name}.finished.json').read_text())
        assert receipt['waited_reaped'] and receipt['returncode']==0
        assert receipt['tsv_sha256']==sha256((D/data).read_bytes()).hexdigest()
        assert receipt['resource_sha256']==sha256((D/f'{name}.resource').read_bytes()).hexdigest()
        assert receipt['source_sha256']==sha256((D/'contact_geometry.cpp').read_bytes()).hexdigest()
    independent,gate=brute(3)
    assert independent==actual[3]
    gate['all_path_weight_cases']=all_path_bound()
    print(json.dumps({'gate':gate,'sizes':[summary(n,actual[n],old[n]) for n in (3,4,5)]},indent=2))
