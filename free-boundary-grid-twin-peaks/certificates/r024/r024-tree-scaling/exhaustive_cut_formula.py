#!/usr/bin/env python3
"""All labelled trees through six vertices: cut formula versus subset DP."""
from itertools import product
from math import factorial,comb
import importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    assert s.loader is not None;s.loader.exec_module(m);return m
t=load('frozen',HERE/'tree_peak_dp_r024_frozen.py')
e=load('subset','/Volumes/data/CodexTemp/T2/a022-global-top2-proof-r024/exact_kernel.py')
def prufer(code):
    n=len(code)+2;deg=[1]*n
    for x in code:deg[x]+=1
    a=[set() for _ in range(n)]
    for x in code:
        u=next(i for i,d in enumerate(deg) if d==1);a[u].add(x);a[x].add(u);deg[u]-=1;deg[x]-=1
    u=[i for i,d in enumerate(deg) if d==1];a[u[0]].add(u[1]);a[u[1]].add(u[0])
    return a
def rooted(a,V,r):
    def size(v,p):return 1+sum(size(w,v) for w in a[v] if w!=p and w in V)
    z=1
    def go(v,p):
        nonlocal z
        for w in a[v]:
            if w!=p and w in V:go(w,v);z*=size(w,v)
    go(r,-1);return factorial(len(V)-1)//z
def cuts(a):
    N=len(a);z=0
    for u in range(N):
        for v in a[u]:
            if u>=v:continue
            A={u};q=[u]
            while q:
                x=q.pop()
                for w in a[x]:
                    if {x,w}=={u,v}:continue
                    if w not in A:A.add(w);q.append(w)
            C=set(range(N))-A
            z+=comb(N-2,len(A)-1)*(sum(rooted(a,A,x) for x in A)*sum(rooted(a,C,x) for x in C)-rooted(a,A,u)*rooted(a,C,v))
    return z
for n in range(3,7):
    done=0
    for code in product(range(n),repeat=n-2):
        a=prufer(code);_,B=e.exact_counts([sum(1<<w for w in row) for row in a]);assert B==cuts(a)
        done+=1
    print(n,done,'passed')
