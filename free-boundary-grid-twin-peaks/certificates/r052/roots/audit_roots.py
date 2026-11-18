"""Check the prescribed-root inverse on n=3 and independently generated words."""
from collections import Counter
from itertools import permutations
from math import perm
from random import Random


def grid(n):
    points=list(range(n*n))
    return [{v for v in points if abs(u//n-v//n)+abs(u%n-v%n)==1} for u in points]


def peaks(word, adj):
    earlier=set()
    result=[]
    for v in word:
        if not earlier & adj[v]:
            result.append(v)
        earlier.add(v)
    return tuple(result)


def patch(n, adj, p, q, x, y):
    def distance(u,v):
        return abs(u//n-v//n)+abs(u%n-v%n)
    vertices={x,y}
    for u in (p,q):
        center=min((x,y),key=lambda v:(distance(u,v),v))
        while u!=center:
            vertices.add(u)
            u=min(v for v in adj[u] if distance(v,center)<distance(u,center))
    ordered=[x,y]
    left=vertices-{x,y}
    seen={x,y}
    while left:
        u=min(v for v in left if adj[v] & seen)
        ordered.append(u)
        left.remove(u)
        seen.add(u)
    return tuple(ordered)


def encode(word, movable, p, adj, x, y):
    P=set(movable)
    rank={v:i for i,v in enumerate(word)}
    carried={rank[v]:v for v in P}
    assert carried[0]==p
    output=movable+tuple(v for v in word if v not in P)
    assert peaks(output,adj)==(x,y)
    other=iter(output[len(movable):])
    reassembled=tuple(carried[i] if i in carried else next(other) for i in range(len(word)))
    assert reassembled==word
    return output


def exhaustive():
    n=3
    target=(0,8)
    adj=grid(n)
    cache={(p,q):patch(n,adj,p,q,*target) for p in range(9) for q in range(9) if p!=q and q not in adj[p]}
    preimages=Counter()
    target_count=0
    for word in permutations(range(9)):
        roots=peaks(word,adj)
        if len(roots)!=2:
            continue
        p,q=roots
        target_count+=set(roots)==set(target)
        output=encode(word,cache[p,q],p,adj,*target)
        preimages[output]+=1
    cap=sum(perm(8,len(P)-1) for P in cache.values())
    assert sum(preimages.values())==143808
    assert 0<len(preimages)<=target_count and max(preimages.values())<=cap
    print('n=3',sum(preimages.values()),target_count,len(preimages),max(preimages.values()),cap)


def samples():
    rng=Random(2026111852)
    total=0
    for n in (4,5,8,12):
        adj=grid(n)
        for _ in range(1000):
            def choose_pair():
                p=rng.randrange(n*n)
                q=rng.choice(sorted(set(range(n*n))-adj[p]-{p}))
                return p,q
            p,q=choose_pair()
            seen={p,q}
            order=[p,q]
            frontera=(adj[p]|adj[q])-seen
            while frontera:
                u=rng.choice(sorted(frontera))
                order.append(u)
                seen.add(u)
                frontera.remove(u)
                frontera |= adj[u]-seen
            assert peaks(order,adj)==(p,q)
            x,y=choose_pair()
            movable=patch(n,adj,p,q,x,y)
            encode(tuple(order),movable,p,adj,x,y)
            assert len(patch(n,adj,p,q,0,n*n-1))<=2*n
            total+=1
    print('seed',2026111852,'constructive_word_target_pairs',total)


if __name__=='__main__':
    exhaustive()
    samples()
