#!/usr/bin/env python3
"""Exact degree-nine support sieve, root envelopes, and singleton trace boxes."""

import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
import json
from math import isqrt, prod
from pathlib import Path

from singleton_probe import house


DEGREE=9


def supports():
    def reason(r,t):
        if r==t:return "equal-first-offset"
        if min(r,t)==1 and max(r,t) in (2,3,4):return "one-versus-two-through-four"
        if {r,t}=={2,4}:return "two-versus-four"
        return None

    allowed=[]
    rejected=Counter()
    for mask in range(1,1<<(DEGREE-1)):
        S=tuple(j for j in range(1,DEGREE) if mask>>(j-1)&1)
        code=reason(S[0],DEGREE-S[-1])
        if code:rejected[code]+=1
        else:allowed.append(S)
    # Independent endpoint expansion: choose the endpoints, then any interior subset.
    banned={(1,2),(1,3),(1,4),(2,4)}
    endpoint_expansion=set()
    for left in range(1,DEGREE):
        for right in range(left,DEGREE):
            pair=tuple(sorted((left,DEGREE-right)))
            if pair[0]==pair[1] or pair in banned:continue
            if left==right:endpoint_expansion.add((left,))
            else:
                interior=range(left+1,right)
                for k in range(len(interior)+1):
                    for inside in combinations(interior,k):
                        endpoint_expansion.add((left,*inside,right))
    assert set(allowed)==endpoint_expansion
    def reversal(S):
        return tuple(DEGREE-j for j in reversed(S))
    representatives=sorted((S for S in allowed if S<reversal(S)),key=lambda s:(len(s),s))
    assert not any(S==reversal(S) for S in allowed)
    assert len(allowed)==2*len(representatives)==50
    return rejected,representatives,sorted(allowed)


def envelopes():
    vertices=set()
    for free in range(DEGREE):
        for sides in product((-1,1),repeat=DEGREE-1):
            last=-sum(sides)
            if abs(last)<=1:
                vertices.add((*sides[:free],last,*sides[free:]))
    assert len(vertices)==630
    assert {tuple(sorted(v)) for v in vertices}=={(-1,-1,-1,-1,0,1,1,1,1)}
    radii=(Q(1,2),)*4+(Q(1),)+(Q(2),)*4
    e=[sum(prod(side) for side in combinations(radii,j)) for j in range(DEGREE+1)]
    T=[sum(x**j for x in radii) for j in range(1,DEGREE+1)]
    b=[Q(1)]
    for x in radii:
        updated=[Q(0)]*(len(b)+1)
        for k,c in enumerate(b):
            updated[k]+=c
            updated[k+1]+=x*c
        b=updated
    assert e==b and e==e[::-1]
    assert T==[4*(Q(2)**j+Q(1,2)**j)+1 for j in range(1,DEGREE+1)]
    return len(vertices),e,T


def domain_box(r,t,T):
    gap=2*min(T[r-1]/r,T[t-1]/t)
    square=gap*gap
    return [D for D in range(2,(square.numerator-1)//square.denominator+1)
            if (1 if D%4==1 else 4)*D<square
            and all(D%(a*a) for a in range(2,isqrt(D)+1))]


def delta_values(T,center,radius):
    margin=radius-abs(T-center)
    if margin<=0:return
    square=margin*margin
    for delta in range(T*T%4,(square.numerator-1)//square.denominator+1,4):
        if delta and isqrt(delta)**2!=delta:
            yield delta


def cap(bound):
    return (bound.numerator-1)//bound.denominator


def first_boxes(e,M):
    r1, r2, r4 = Counter(), Counter(), Counter()
    Lbound=cap(e[1])
    p2bound=cap(e[2])
    for T in range(-cap(2*e[1]),cap(2*e[1])+1):
        for d in delta_values(T,0,2*e[1]):
            r1["trace_norm"]+=1
            N=(T*T-d)//4
            p2s=[c for c in range(-p2bound,p2bound+1)
                 if house((T*T-2*N-4*c,T),M[1].numerator,M[1].denominator,d)]
            if p2s:
                r1["trace_norm_with_newton_2"]+=1
                r1["first_two_coefficients"]+=len(p2s)
    for L in range(-Lbound,Lbound+1):
        centre=L*L
        radius=M[1]
        for T in range(centre-cap(radius),centre+cap(radius)+1):
            for d in delta_values(T,centre,radius):
                if house((T,1),e[2].numerator,e[2].denominator,d):
                    r2["first_two_coefficients"]+=1
    p2max=cap(M[1]/2)
    p3max=cap(M[2]/3)
    radius=M[3]/2
    for p2 in range(-p2max,p2max+1):
        centre=p2*p2
        for T in range(centre-cap(radius),centre+cap(radius)+1):
            for d in delta_values(T,centre,radius):
                if house((T,1),e[4].numerator,e[4].denominator,d):
                    r4["p2_and_trace_norm"]+=1
                    r4["first_four_coefficients"]+=2*p3max+1
    return {1:r1,2:r2,4:r4}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert not a.out.exists()
    forbidden,reps,all_live=supports()
    n,e,T=envelopes()
    pairs=sorted({(S[0],DEGREE-S[-1]) for S in reps})
    fields={f"{r},{t}": domain_box(r,t,T) for r,t in pairs}
    report=dict(status="verified",vertices=n,coefficient_bound=list(map(str,e)),
                power_sum_bound=list(map(str,T)),supports_examined=255,
                endpoint_rejections=forbidden,ordinary_supports=all_live,
                reversal_representatives=reps,field_counts={k:len(v) for k,v in fields.items()},
                field_lists=fields,single_site_early_boxes=first_boxes(e,T))
    a.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k: report[k] for k in ("vertices","coefficient_bound","power_sum_bound",
                     "endpoint_rejections","field_counts","single_site_early_boxes")},sort_keys=True))


if __name__=="__main__":
    main()
