#!/usr/bin/env python3
"""Exact n=3 and seeded light n=4 checks of the reinsertion harmonic sum."""
from fractions import Fraction
from itertools import permutations
import json, random
from pathlib import Path


def grid(n):
    N=n*n; adj=[set() for _ in range(N)]
    for r in range(n):
        for c in range(n):
            v=r*n+c
            for dr,dc in ((1,0),(-1,0),(0,1),(0,-1)):
                rr,cc=r+dr,c+dc
                if 0<=rr<n and 0<=cc<n: adj[v].add(rr*n+cc)
    return adj


def is_onepeak(word,adj):
    seen=set()
    for i,v in enumerate(word):
        if i and not (adj[v]&seen): return False
        seen.add(v)
    return True


def harmonic(word,adj):
    N=len(word); pos=[0]*N
    for i,v in enumerate(word): pos[v]=i
    r=word[0]; out=Fraction(0)
    rs=[]
    for p in range(N):
        if p==r or p in adj[r]: continue
        lp=min(pos[u] for u in adj[p])
        bp=N
        for rank in range(pos[p]+1,N):
            v=word[rank]
            early=[u for u in adj[v] if pos[u]<rank]
            if early==[p]:
                bp=rank; break
        # zero-based indices: valid insertion positions are lp+1,...,bp-1,
        # so R=bp-lp-1, equal to 1-based b_p-l_p-1.
        R=bp-lp-1
        assert R>=1
        rs.append([p,R])
        out += Fraction(1,R)
    return out,rs


def n3_exact():
    adj=grid(3); N=9; num=0; best=None; witnesses=[]; total=Fraction(0)
    for w in permutations(range(N)):
        if not is_onepeak(w,adj): continue
        num+=1; val,rs=harmonic(w,adj); total+=val
        if best is None or val<best:
            best=val; witnesses=[[list(w),rs]]
        elif val==best and len(witnesses)<5:
            witnesses.append([list(w),rs])
    return {"complete":True,"onepeak_words":num,"minimum":[best.numerator,best.denominator],
            "sum_harmonic":[total.numerator,total.denominator],
            "mean_harmonic":[(total/num).numerator,(total/num).denominator],
            "witnesses":witnesses}


def n4_seeded():
    adj=grid(4); N=16; rng=random.Random(20261117); best=None; witness=None
    # 4096 lightweight growth histories, not an exhaustive n=4 enumeration.
    for _ in range(4096):
        r=rng.randrange(N); seen={r}; w=[r]
        while len(w)<N:
            front=sorted(set().union(*(adj[v] for v in seen))-seen)
            v=front[rng.randrange(len(front))]
            seen.add(v); w.append(v)
        assert is_onepeak(w,adj)
        val,rs=harmonic(w,adj)
        if best is None or val<best:
            best=val; witness=[w,rs]
    return {"complete":False,"method":"seeded 4096 growth histories; no minimum certificate",
            "seed":20261117,"minimum_seen":[best.numerator,best.denominator],"witness":witness}


out={"definition":"R=b_p-l_p-1 in one-based ranks; code uses equal zero-based R=bp-lp-1",
     "n3":n3_exact(),"n4":n4_seeded()}
Path(__file__).with_name("reinsertion-sum-small.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(json.dumps(out,indent=2,sort_keys=True))
