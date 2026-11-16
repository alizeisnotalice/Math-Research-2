#!/usr/bin/env python3
"""Frozen r024 exact low-peak and top-two-root probabilities on trees."""
from fractions import Fraction as F
from math import factorial
ZMAX=2
def trim(p):
    while p and not p[-1]:p.pop()
    return p
def add(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,x in enumerate(p):r[i]+=x
    for i,x in enumerate(q):r[i]+=x
    return trim(r)
def sub(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,x in enumerate(p):r[i]+=x
    for i,x in enumerate(q):r[i]-=x
    return trim(r)
def mul(p,q):
    if not p or not q:return []
    r=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):r[i+j]+=x*y
    return trim(r)
def integ0(p):return [F(0)]+[x/F(i+1) for i,x in enumerate(p)]
def integ1(p):return sum((x/F(i+1) for i,x in enumerate(p)),F(0))
def const(c):return [F(c)] if c else []
ONE=const(1)
def zp_zero():return [[] for _ in range(ZMAX+1)]
def zp_one():return [ONE]+[[] for _ in range(ZMAX)]
def zp_add(a,b):return [add(a[i],b[i]) for i in range(ZMAX+1)]
def zp_mul(a,b):
    out=zp_zero()
    for i in range(ZMAX+1):
        for j in range(ZMAX+1-i):out[i+j]=add(out[i+j],mul(a[i],b[j]))
    return out
def zp_int0(a):return [integ0(x) for x in a]
def zp_int1(a):return [integ1(x) for x in a]
def zp_const(v):return [const(x) for x in v]
def zp_zminus1(a):
    out=zp_zero()
    for i in range(ZMAX+1):
        out[i]=sub(out[i],a[i])
        if i+1<=ZMAX:out[i+1]=add(out[i+1],a[i])
    return out
def zp_coef_integral(a,k):return integ1(a[k])
def peak_probability_poly(adj,root=0):
    def rec(v,p):
        A=zp_one();B=zp_one()
        for w in adj[v]:
            if w==p:continue
            Aw,Bw=rec(w,v);tail=[sub(const(integ1(poly)),integ0(poly)) for poly in Bw]
            A=zp_mul(A,zp_add(zp_const(zp_int1(Aw)),zp_zminus1(tail)))
            B=zp_mul(B,zp_int0(Aw))
        return A,B
    A,B=rec(root,-1);final=zp_add(A,zp_zminus1(B))
    return [zp_coef_integral(final,k) for k in range(ZMAX+1)]
def no_peak_component(adj,vertices,marked):
    start=next(v for v in vertices if v in marked);seen=set()
    def rec(v,p):
        seen.add(v);A=ONE;B=ONE
        for w in adj[v]:
            if w==p or w not in vertices:continue
            Aw,Bw=rec(w,v)
            child=const(integ1(Aw)) if w in marked else sub(const(integ1(Aw)),sub(const(integ1(Bw)),integ0(Bw)))
            A=mul(A,child);B=mul(B,integ0(Aw))
        return A,B
    A,_=rec(start,-1);assert seen==vertices;return integ1(A)
def top2_probability(adj):
    n=len(adj);total=F(0)
    for x in range(n):
        for y in range(n):
            if x==y or y in adj[x]:continue
            rem=set(range(n))-{x,y};q=F(1)
            while rem:
                seed=next(iter(rem));stack=[seed];comp={seed};rem.remove(seed)
                while stack:
                    v=stack.pop()
                    for w in adj[v]:
                        if w in rem:rem.remove(w);comp.add(w);stack.append(w)
                marked={v for v in comp if x in adj[v] or y in adj[v]};assert marked
                q*=no_peak_component(adj,comp,marked)
            total+=q
    return total/F(n*(n-1))
def binary(h):
    n=(1<<(h+1))-1;a=[set() for _ in range(n)]
    for i in range(1,n):a[i].add((i-1)//2);a[(i-1)//2].add(i)
    return a
def double_binary(h):
    a=binary(h);m=len(a);g=[set(x) for x in a]+[{m+w for w in x} for x in a]
    g[0].add(m);g[m].add(0);return g
def binary_orbit_groups(adj):
    n=len(adj);dep=[(v+1).bit_length()-1 for v in range(n)];par=[(v-1)//2 if v else -1 for v in range(n)]
    def lca(a,b):
        while dep[a]>dep[b]:a=par[a]
        while dep[b]>dep[a]:b=par[b]
        while a!=b:a=par[a];b=par[b]
        return a
    groups={}
    for x in range(n):
        for y in range(n):
            if x==y or y in adj[x]:continue
            groups.setdefault((dep[x],dep[y],dep[lca(x,y)]),[0,(x,y)])[0]+=1
    return groups
def pairq(adj,x,y):
    rem=set(range(len(adj)))-{x,y};q=F(1)
    while rem:
        seed=next(iter(rem));stack=[seed];comp={seed};rem.remove(seed)
        while stack:
            v=stack.pop()
            for w in adj[v]:
                if w in rem:rem.remove(w);comp.add(w);stack.append(w)
        marked={v for v in comp if x in adj[v] or y in adj[v]};assert marked
        q*=no_peak_component(adj,comp,marked)
    return q
def top2_probability_binary_fast(adj):
    n=len(adj);groups=binary_orbit_groups(adj)
    return sum((count*pairq(adj,*rep) for count,rep in groups.values()),F(0))/F(n*(n-1)),len(groups)
def double_binary_orbit_contributions(h):
    adj=double_binary(h);m=(1<<(h+1))-1;n=2*m;dep=[(v+1).bit_length()-1 for v in range(m)];par=[(v-1)//2 if v else -1 for v in range(m)]
    def lca(a,b):
        while dep[a]>dep[b]:a=par[a]
        while dep[b]>dep[a]:b=par[b]
        while a!=b:a=par[a];b=par[b]
        return a
    groups={}
    for x in range(n):
        for y in range(x+1,n):
            if y in adj[x]:continue
            bx,by=x//m,y//m;u=x%m;v=y%m
            key=('S',min(dep[u],dep[v]),max(dep[u],dep[v]),dep[lca(u,v)]) if bx==by else ('X',min(dep[u],dep[v]),max(dep[u],dep[v]))
            groups.setdefault(key,[0,(x,y)])[0]+=1
    return adj,[(key,count,rep,pairq(adj,*rep)) for key,(count,rep) in groups.items()]
def top2_probability_double_binary_fast(h):
    adj,rows=double_binary_orbit_contributions(h);n=len(adj)
    return 2*sum((count*q for _,count,_,q in rows),F(0))/F(n*(n-1)),len(rows)
def counts(adj):
    n=len(adj);p=peak_probability_poly(adj);a2=p[2]*factorial(n);b=top2_probability(adj)*factorial(n)
    assert a2.denominator==b.denominator==1;return int(a2),int(b),p
def binary_counts_fast(h):
    adj=binary(h);n=len(adj);a2=peak_probability_poly(adj)[2]*factorial(n);b,g=top2_probability_binary_fast(adj);b*=factorial(n)
    assert a2.denominator==b.denominator==1;return int(a2),int(b),g
def double_binary_counts_fast(h):
    adj=double_binary(h);n=len(adj);a2=peak_probability_poly(adj)[2]*factorial(n);b,g=top2_probability_double_binary_fast(h);b*=factorial(n)
    assert a2.denominator==b.denominator==1;return int(a2),int(b),g
def path(n):
    a=[set() for _ in range(n)]
    for i in range(n-1):a[i].add(i+1);a[i+1].add(i)
    return a
def rooted_onepeak_count(adj):
    n=len(adj);ans=0
    for r in range(n):
        def size(v,p):return 1+sum(size(w,v) for w in adj[v] if w!=p)
        prod=1
        def walk(v,p):
            nonlocal prod
            for w in adj[v]:
                if w!=p:walk(w,v)
            if v!=r:prod*=size(v,p)
        walk(r,-1);ans+=factorial(n-1)//prod
    return ans
