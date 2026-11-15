"""Complete 3x3 audit of the internal-shortest-path lifting injection."""
from collections import Counter
from itertools import permutations
from math import prod
import json
from exact import grid,peaks


def path(n,u,v):
    x,y=u%n,u//n;X,Y=v%n,v//n;out=[u]
    while x!=X:x+=1 if x<X else -1;out.append(y*n+x)
    while y!=Y:y+=1 if y<Y else -1;out.append(y*n+x)
    return out


def main():
    n=3;N=n*n;_,adj=grid(n,n);fibers=Counter();distance=Counter()
    for L in permutations(range(1,N+1)):
        pp=peaks(L,adj)
        if len(pp)!=2:continue
        u=max(pp,key=lambda v:L[v]);v=min(pp,key=lambda v:L[v]);P=path(n,u,v);d=len(P)-1
        inside=P[1:-1]
        order=[u]+inside+[w for w in sorted(range(N),key=lambda w:-L[w]) if w!=u and w not in inside]
        out=[0]*N
        for rank,w in enumerate(order):out[w]=N-rank
        assert peaks(out,adj)==[u]
        assert order[:d]==P[:-1]
        # Restore the old internal ranks; the complement has its old relative order.
        slots={w:L[w] for w in inside};remaining=sorted(set(range(1,N+1))-set(slots.values()),reverse=True)
        restored=[0]*N
        for w,value in slots.items():restored[w]=value
        for w,value in zip([w for w in order if w not in slots],remaining):restored[w]=value
        assert tuple(restored)==L
        fibers[d,tuple(out)]+=1;distance[d]+=1
    maximum={d:max(c for (e,_),c in fibers.items() if e==d) for d in distance}
    assert all(maximum[d]<=prod(range(N-d+1,N)) for d in distance)
    print(json.dumps({'n':n,'checked_A2':sum(distance.values()),'distance_counts':dict(distance),
                      'max_fiber_by_distance':maximum,'universal_fiber_bounds':{d:prod(range(N-d+1,N)) for d in distance}},indent=2))


if __name__=='__main__':main()
