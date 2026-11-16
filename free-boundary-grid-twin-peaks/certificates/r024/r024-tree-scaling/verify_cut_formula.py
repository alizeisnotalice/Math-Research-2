#!/usr/bin/env python3
"""Independent hook-length verification of the tree cut formula for B."""
from math import factorial,comb
from fractions import Fraction
import importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen',HERE/'tree_peak_dp_r024_frozen.py')
t=importlib.util.module_from_spec(spec);assert spec.loader is not None;spec.loader.exec_module(t)

def rooted_count(adj,vertices,root):
    def size(v,p):return 1+sum(size(w,v) for w in adj[v] if w!=p and w in vertices)
    product=1
    def walk(v,p):
        nonlocal product
        for w in adj[v]:
            if w!=p and w in vertices:walk(w,v);product*=size(w,v)
    walk(root,-1)
    return factorial(len(vertices)-1)//product

def B_by_cuts(adj):
    N=len(adj);total=0
    for u in range(N):
        for v in adj[u]:
            if u>=v:continue
            A={u};todo=[u]
            while todo:
                x=todo.pop()
                for y in adj[x]:
                    if {x,y}=={u,v}:continue
                    if y not in A:A.add(y);todo.append(y)
            C=set(range(N))-A
            A1=sum(rooted_count(adj,A,r) for r in A)
            C1=sum(rooted_count(adj,C,r) for r in C)
            total+=comb(N-2,len(A)-1)*(A1*C1-rooted_count(adj,A,u)*rooted_count(adj,C,v))
    return total

def hook_split_factor_check(adj,r=0):
    """Verify corrected (C2), including ancestor hook distortion G."""
    N=len(adj)
    def root_parent_sizes(root):
        par={root:-1};order=[root]
        for x in order:
            for y in adj[x]:
                if y!=par[x]:par[y]=x;order.append(y)
        sz={v:1 for v in order}
        for v in reversed(order[1:]):sz[par[v]]+=sz[v]
        return par,sz
    par,sz=root_parent_sizes(r)
    ar=rooted_count(adj,set(range(N)),r)
    for u in range(N):
        if u==r:continue
        v=par[u]
        A={u};todo=[u]
        while todo:
            x=todo.pop()
            for y in adj[x]:
                if {x,y}=={u,v}:continue
                if y not in A:A.add(y);todo.append(y)
        C=set(range(N))-A;s=len(A)
        G=1
        w=v
        while w!=r:
            G*=Fraction(sz[w]-s,sz[w]);w=par[w]
        rhs=Fraction(comb(N-1,s)*rooted_count(adj,A,u)*rooted_count(adj,C,r))*G
        assert rhs==ar,(r,u,v,ar,rhs,G)
    return ar

cases={'P4':t.path(4),'P5':t.path(5),'binary-h3':t.binary(3),
       'double-h2':t.double_binary(2),'double-h3':t.double_binary(3),
       'double-h4':t.double_binary(4)}
for name,g in cases.items():
    _,b,_=t.counts(g);c=B_by_cuts(g);assert b==c,(name,b,c)
    print(name,'B=',b,'hook-split=',hook_split_factor_check(g),'passed')
