#!/usr/bin/env python3
"""Exact n=3,5 audit of the complete-square two-sided completion obstruction."""
from itertools import permutations
from math import factorial
import json
from pathlib import Path


def make(n):
    m=(n-1)//2
    V={(x,y) for x in range(n) for y in range(n)}
    A={(x,y) for x in range(m) for y in range(n-1)}
    B={(x,y) for x in range(m+1,n) for y in range(n-1)}
    C=V-A-B
    adj={v:set() for v in V}
    for x,y in V:
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            q=x+dx,y+dy
            if q in V:adj[x,y].add(q)
    return m,V,A,B,C,adj


def count(seed, remaining, adj):
    z=0
    for word in permutations(remaining):
        seen=set(seed)
        for v in word:
            if not adj[v]&seen:break
            seen.add(v)
        else:z+=1
    return z


def connected(seed, vertices, adj):
    seen={seed};todo=[seed]
    while todo:
        v=todo.pop()
        for q in adj[v]&vertices:
            if q not in seen:seen.add(q);todo.append(q)
    return seen==vertices


rows=[]
for n in (3,5):
    m,V,A,B,C,adj=make(n);g=len(C)
    assert connected(next(iter(A)),A,adj) and connected(next(iter(B)),B,adj)
    assert connected(next(iter(B)),V-A,adj) and connected(next(iter(A)),V-B,adj)
    assert not any(adj[x]&B for x in A)
    W=[w for w in C if adj[w]&A and adj[w]&B]
    ha=count(A,tuple(C),adj);hb=count(B,tuple(C),adj)
    f=sum(count(A|B|{w},tuple(C-{w}),adj) for w in W)
    special=(m,n-2)
    special_h=count(A|B|{special},tuple(C-{special}),adj)
    upper=factorial(g)//factorial(m+1)
    assert special in W and special_h==factorial(g-1)
    assert ha<=upper and hb<=upper
    rows.append({'n':n,'m':m,'C_size':g,'common_bridges':W,
                 'H_A':ha,'H_B':hb,'H_upper':upper,'F_all_bridges':f,
                 'special_bridge':special,'special_completion':special_h,
                 'lower_ratio_F_over_maxH':f'{special_h}/{max(ha,hb)}',
                 'symbolic_lower_ratio':f'{factorial(m+1)}/{g}'})
out={'status':'passed','rows':rows}
Path('free-square-obstruction.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
