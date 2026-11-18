#!/usr/bin/env python3
"""Complete nine-degree r=3 singleton mother entry and exact root screening."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
from functools import lru_cache
from math import isqrt
from pathlib import Path

from r3_entry import core_states
from singleton_probe import (ONE, ZERO, ci, conj, plus, scale, house,
                             newton_base, prefixes, root_code, extend, coeff_base)

BOUND=[None]+[4*(Fraction(2)**k+Fraction(1,2)**k)+1 for k in range(1,10)]
COEFF=[1,10,51,133,212,212,133,51,10,1]


def shift_set(base,k,D):
    u,v=base
    n,d=BOUND[k].numerator,BOUND[k].denominator
    margin=(2*n-isqrt(D*(d*v)**2)-1)//d
    if margin<0:return ()
    return range(max(-COEFF[k],-((margin-u)//(2*k))),min(COEFF[k],(margin+u)//(2*k))+1)


def passes(s,k,D):
    return house(s,BOUND[k].numerator,BOUND[k].denominator,D)


@lru_cache
def reverse_pair(e):
    return [(g,h) for h in range(-10,11) for g in range(-51,52)
            if abs(h*h-2*g)<18 and 2*abs(-h**3+3*h*g-3*e)<67]


def fraction_rows(A,B,a,G,Q,k,D):
    return [(A+(x,),B+(y,)) for _,_,x,y in extend(A,B,scale(G,a),scale(Q,conj(a)),k,D)]


def changes(G,L,M,c,d,e):
    g4=-c
    g5=-d-c*G[1]-L*g4
    g6=-e-c*G[2]-d*G[1]-L*g5-M*g4
    return G[4]+g4,G[5]+g5,G[6]+g6


def reverse_compatible(R,Lr,Ls,D):
    Rb=tuple(map(conj,R))
    states=[((1,),(1,),(ONE,),(ONE,))]
    for k in range(1,10):
        nxt=[]
        for U,V,A,B in states:
            for u,v,a,b in extend(A,B,coeff_base(R,U,k),coeff_base(Rb,V,k),k,D):
                if k==3 and (a!=ZERO or b!=ZERO):continue
                if (a!=ZERO and Ls[k]!=ZERO) or (b!=ZERO and Lr[k]!=ZERO):continue
                nxt.append((U+(u,),V+(v,),A+(a,),B+(b,)))
        if not nxt:return False
        states=nxt
    return True


def solve(core):
    L,M,T,D=core
    a=(T,1)
    P3=(ONE,ci(L),ci(M),a)
    S3=[ci(9)]
    for k in range(1,4):
        S3.append(plus(newton_base(P3,S3,k,D),scale(-k,P3[k])))
    low=prefixes(P3+(ZERO,)*3,D,3,6,True)
    counts=Counter()
    if not low:return counts,[]
    assert len(low)==1
    G,Q,A,B=low[0]
    counts["row6_cores"]+=1
    rows=[]
    for c in shift_set(newton_base(P3,S3,4,D),4,D):
        g4,q4=G[4]-c,Q[4]-c
        states7=fraction_rows(A,B,a,g4,q4,7,D)
        if not states7:continue
        counts["row7_states"]+=len(states7)
        P4=P3+(ci(c),)
        S4=S3+[plus(newton_base(P3,S3,4,D),ci(-4*c))]
        base5=newton_base(P4,S4,5,D)
        for d in shift_set(base5,5,D):
            _,g5,_=changes(G,L,M,c,d,0)
            _,q5,_=changes(Q,L,M,c,d,0)
            states8=[out for X,Y in states7 for out in fraction_rows(X,Y,a,g5,q5,8,D)]
            if not states8:continue
            counts["row8_states"]+=len(states8)
            P5=P4+(ci(d),)
            S5=S4+[plus(base5,ci(-5*d))]
            base6=newton_base(P5,S5,6,D)
            for e in shift_set(base6,6,D):
                _,_,g6=changes(G,L,M,c,d,e)
                _,_,q6=changes(Q,L,M,c,d,e)
                states9=[out for X,Y in states8 for out in fraction_rows(X,Y,a,g6,q6,9,D)]
                if not states9:continue
                counts["row9_states"]+=len(states9)
                P6=P5+(ci(e),)
                S6=S5+[plus(base6,ci(-6*e))]
                base7=newton_base(P6,S6,7,D)
                gs=shift_set(base7,7,D)
                for g,h in reverse_pair(e):
                    if g not in gs:continue
                    P=P6+(ci(g),ci(h),ONE)
                    s7=plus(base7,ci(-7*g))
                    s8=plus(newton_base(P,S6+[s7],8,D),ci(-8*h))
                    if not passes(s8,8,D):continue
                    s9=plus(newton_base(P,S6+[s7,s8],9,D),ci(-9))
                    if not passes(s9,9,D):continue
                    R=P[::-1]
                    S=[ci(9)]
                    for k in range(1,10):
                        sk=plus(newton_base(R,S,k,D),scale(-k,R[k]))
                        if not passes(sk,k,D):break
                        S.append(sk)
                    if len(S)!=10:continue
                    counts["newton_both"]+=1
                    if not any(reverse_compatible(R,X,Y,D) for X,Y in states9):continue
                    counts["coupled_mothers"]+=1
                    code=root_code(P,D)
                    rows.append(((L,M,T,D,c,d,e,g,h),code))
                    counts[code]+=1
    return counts,rows


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    assert not a.out.exists()
    cores=[row[:4] for row in core_states()]
    totals=Counter()
    rows=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for stats,records in pool.map(solve,cores,chunksize=1):
            totals.update(stats)
            rows+=records
    assert len(rows)==len(set(row for row,_ in rows))
    a.out.with_suffix(".rows.jsonl").write_text(
        "".join(json.dumps(x,separators=(",",":"))+"\n" for x in sorted(rows)))
    result=dict(status="complete_entry",input_cores=len(cores),counts=totals)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
