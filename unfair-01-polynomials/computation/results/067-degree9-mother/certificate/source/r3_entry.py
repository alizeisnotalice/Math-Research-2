#!/usr/bin/env python3
"""Exhaust the first fractional core and the row-six budget for the r=3 site."""

import argparse
from collections import Counter
from itertools import product
import json
from math import isqrt
from pathlib import Path

from singleton_probe import ONE, ZERO, ci, house, prefixes


def core_states():
    for L in range(-10,11):
        for M in range(-51,52):
            if abs(L*L-2*M)>=18:continue
            base=-L**3+3*L*M
            allowed={}
            for (x,y),(a,b) in product(((0,0),(1,0),(0,1)),repeat=2):
                if x*b or y*a:continue
                w=L**3-2*L*M+L*a+(M-L*L)*x
                z=L**3-2*L*M+L*b+(M-L*L)*y
                T=1-w-z
                margin=67-abs(2*base-3*T)
                if margin<=0:continue
                upper=(margin*margin-1)//9
                for disc in range(T*T%4,upper+1,4):
                    if not disc or isqrt(disc)**2==disc:continue
                    if not house((T,1),134,1,disc):continue
                    shift=T+2*w
                    if shift>=2:continue
                    if shift<0 and disc<=shift*shift:continue
                    if disc>=(2-shift)**2:continue
                    allowed.setdefault((T,disc),set()).add((x,y,a,b))
            for (T,disc),low in sorted(allowed.items()):
                yield L,M,T,disc,low


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert not a.out.exists()
    counts=Counter()
    rows=[]
    for L,M,T,disc,low in core_states():
        counts["distinct_first_fraction_cores"]+=1
        counts["row3_states"]+=len(low)
        P=(ONE,ci(L),ci(M),(T,1),ZERO,ZERO,ZERO)
        direct=prefixes(P,disc,3,6,True)
        if direct:
            counts["distinct_row6_cores"]+=1
            counts["row6_states"]+=len(direct)
            rows.append((L,M,T,disc,len(low),len(direct)))
    report=dict(counts=counts,rows=rows,
                scope="Complete Newton and first fractional cores for r=3; all low row-six budgets after absorbing rational p4,p5,p6 into integer shifts.")
    a.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps(counts,sort_keys=True))


if __name__=="__main__":
    main()
