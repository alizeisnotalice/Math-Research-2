#!/usr/bin/env python3
"""Exhaustive n4 verification of the corrected fibre bound b-3."""
from itertools import permutations
from collections import Counter
import json
n=4
V={(x,y) for x in range(n) for y in range(n)}
S={(2,0),(0,1),(1,1),(2,1),(3,0),(0,2),(0,3)}
A={(0,0),(1,0)}
B=V-S-A
r=(1,2)
def adj(v):
    x,y=v
    return {(x+dx,y+dy) for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]}&V
F=set().union(*(adj(v) for v in S))&B
assert adj(r)&B <=F
core=B-{r}; initial=F-{r}; fixed={x:min(adj(x)&(core-initial)) for x in initial}
def legal(w):
    seen=set()
    for x in w:
        if x not in initial and not adj(x)&seen:return False
        seen.add(x)
    return True
orders=[w for w in permutations(sorted(core)) if legal(w)]
loads=Counter()
for w in orders:
    if w[1] in initial:
        c=fixed[w[0]]
        assert w.index(c)>=2
        image=(w[0],c)+tuple(x for x in w[1:] if x!=c)
        assert legal(image) and image[1] not in initial
        loads[image]+=1
b=len(B);H=len(orders);D=sum(w[1] not in initial for w in orders)
assert H-D<=max(loads.values())*D<=(b-3)*D
assert max(loads.values())<=b-3
print(json.dumps({'n':n,'b':b,'core_orders':H,'core_dependent':D,'bad_orders':H-D,'maximum_fibre':max(loads.values()),'bound':b-3,'status':'passed'},indent=2))
