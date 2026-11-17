#!/usr/bin/env python3
"""Direct exhaustive P3xP3 check of the high/low singleton orientation identity."""
from fractions import Fraction
from itertools import permutations
import hashlib, json, platform
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

def births(word,adj):
    seen=set();out=[]
    for v in word:
        if not(adj[v]&seen):out.append(v)
        seen.add(v)
    return out

def comps(seen,adj):
    todo=set(seen);out=[]
    while todo:
        c={todo.pop()};q=list(c)
        while q:
            u=q.pop(); hit=adj[u]&todo; todo-=hit; c|=hit; q.extend(hit)
        out.append(c)
    return out

def first_merge(word,adj):
    seen=set()
    for v in word:
        cs=comps(seen,adj) if seen else []
        hit=[c for c in cs if c&adj[v]]
        if len(hit)>=2:
            assert len(cs)==len(hit)==2
            return hit
        seen.add(v)
    raise AssertionError

def geom(u,v,n=3):
    dx=abs(u%n-v%n);dy=abs(u//n-v//n)
    return dx+dy,Fraction(dx*dx+dy*dy,n*n)

def main():
    adj=grid(); high_count=0; low_rhs=Fraction()
    high_d=Fraction(); low_d_rhs=Fraction()
    high_E=Fraction(); low_E_rhs=Fraction()
    high_far3=0; low_far3_rhs=Fraction(); k2=0
    for word in permutations(range(9)):
        bs=births(word,adj)
        if len(bs)!=2:continue
        r,p=bs; cs=first_merge(word,adj)
        C_hi=next(c for c in cs if r in c); C_lo=next(c for c in cs if p in c)
        q=len(C_hi)+len(C_lo); d,E=geom(r,p)
        if len(C_hi)==1:
            high_count+=1; high_d+=d; high_E+=E
            if d>=3:high_far3+=1
        if len(C_lo)==1:
            low_rhs+=Fraction(1,q-1)
            low_d_rhs+=Fraction(d,q-1);low_E_rhs+=E/Fraction(q-1)
            if d>=3:low_far3_rhs+=Fraction(1,q-1)
        if q==2:k2+=1
    assert high_count==low_rhs
    assert high_d==low_d_rhs and high_E==low_E_rhs
    assert high_far3==low_far3_rhs
    out={"kind":"direct-P3xP3-high-low-singleton-orientation-audit",
         "python":platform.python_version(),"high_singleton_count":high_count,
         "low_singleton_weighted_count":[low_rhs.numerator,low_rhs.denominator],
         "high_singleton_distance_sum":[high_d.numerator,high_d.denominator],
         "low_singleton_weighted_distance_sum":[low_d_rhs.numerator,low_d_rhs.denominator],
         "high_singleton_E_sum":[high_E.numerator,high_E.denominator],
         "low_singleton_weighted_E_sum":[low_E_rhs.numerator,low_E_rhs.denominator],
         "far_d_ge_3_high_count":high_far3,
         "far_d_ge_3_low_weighted_count":[low_far3_rhs.numerator,low_far3_rhs.denominator],
         "q2_words":k2,"all_identities_passed":True}
    text=json.dumps(out,indent=2,sort_keys=True)+"\n"
    HERE.joinpath("audit-high-singleton-identity-n3.json").write_text(text)
    print(text,end="");print("output_sha256="+hashlib.sha256(text.encode()).hexdigest())
if __name__=="__main__":main()
