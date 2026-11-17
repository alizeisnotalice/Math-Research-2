#!/usr/bin/env python3
"""Exact P3xP3 audit of geometric Psi and singleton endpoint kappa mass."""
from fractions import Fraction
from itertools import permutations
import hashlib,json,platform
from pathlib import Path

HERE=Path(__file__).resolve().parent

def grid(n=3):
    adj=[set() for _ in range(n*n)]
    for y in range(n):
        for x in range(n):
            u=y*n+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                xx,yy=x+dx,y+dy
                if 0<=xx<n and 0<=yy<n:adj[u].add(yy*n+xx)
    return adj

def births(w,adj):
    seen=set();out=[]
    for v in w:
        if not(adj[v]&seen):out.append(v)
        seen.add(v)
    return out

def comps(seen,adj):
    todo=set(seen);out=[]
    while todo:
        c={todo.pop()};q=list(c)
        while q:
            u=q.pop();hit=adj[u]&todo;todo-=hit;c|=hit;q.extend(hit)
        out.append(c)
    return out

def first_merge(w,adj):
    seen=set()
    for v in w:
        cs=comps(seen,adj) if seen else []
        hit=[c for c in cs if c&adj[v]]
        if len(hit)>=2:
            assert len(cs)==len(hit)==2
            return hit
        seen.add(v)
    raise AssertionError

def E(u,v,n=3):
    ux,uy=u%n,u//n;vx,vy=v%n,v//n
    return Fraction((ux-vx)**2+(uy-vy)**2,n*n)

def psi_geo(w,adj):
    N=len(w);pos=[0]*N
    for i,v in enumerate(w,1):pos[v]=i
    r=w[0];total=Fraction()
    for p in range(N):
        if p==r or p in adj[r]:continue
        ell=min(pos[u] for u in adj[p]);b=N+1
        for rank in range(pos[p]+1,N+1):
            v=w[rank-1]
            if [u for u in adj[v] if pos[u]<rank]==[p]:b=rank;break
        R=b-ell-1;assert R>=1
        total+=Fraction(ell-1,R)*E(r,p)
    return total

def main():
    adj=grid(); source=Fraction(); a1=0
    bplus_geo=Fraction();B=0;endpoint=Fraction();exactone=Fraction();both=Fraction()
    for w in permutations(range(9)):
        bs=births(w,adj)
        if len(bs)==1:
            a1+=1;source+=psi_geo(w,adj)
        if len(bs)!=2:continue
        cs=first_merge(w,adj);low=bs[1]
        if len(next(c for c in cs if low in c))==1:
            bplus_geo+=E(bs[0],bs[1])
        if bs!=[w[0],w[1]]:continue
        B+=1
        sizes=sorted(map(len,cs));sing=sum(s==1 for s in sizes)
        if not sing:continue
        a,b=sizes;q=a+b;m=Fraction(q*(q-1),2*a*b);z=E(w[0],w[1])
        endpoint+=m*z
        if sing==2:both+=m*z
        else:exactone+=m*z
    assert a1==29568 and B==81984 and source==bplus_geo
    # On exact-one fibres kappa endpoint mass is q/(q-1) times B+ mass,
    # hence globally it lies between B+ mass and 3B+/2.
    assert bplus_geo<=endpoint<=Fraction(3,2)*bplus_geo
    out={"kind":"exact-P3xP3-geometric-Psi-singleton-endpoint-audit",
         "python":platform.python_version(),"A1":a1,"B":B,
         "sum_geometric_Psi":[source.numerator,source.denominator],
         "direct_Bplus_geometric_mass":[bplus_geo.numerator,bplus_geo.denominator],
         "equality_geometric_Psi_equals_Bplus":source==bplus_geo,
         "raw_B_singleton_endpoint_Ekappa":[endpoint.numerator,endpoint.denominator],
         "raw_exactone_endpoint_Ekappa":[exactone.numerator,exactone.denominator],
         "raw_bothsingleton_endpoint_Ekappa":[both.numerator,both.denominator],
         "endpoint_between_Bplus_and_3over2_Bplus":True}
    text=json.dumps(out,indent=2,sort_keys=True)+'\n';HERE.joinpath('audit-geometric-psi-n3.json').write_text(text)
    print(text,end='');print('output_sha256='+hashlib.sha256(text.encode()).hexdigest())
if __name__=='__main__':main()
