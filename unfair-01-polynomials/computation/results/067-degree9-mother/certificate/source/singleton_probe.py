#!/usr/bin/env python3
"""Exact partial witness search for nine-degree single-site necessary conditions."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path

ZERO, ONE = (0,0), (2,0)
def ci(x):
    return 2*x,0


def plus(x,y):
    return x[0]+y[0],x[1]+y[1]


def neg(x):
    return -x[0],-x[1]


def sub(x,y):
    return x[0]-y[0],x[1]-y[1]


def scale(n,x):
    return n*x[0],n*x[1]


def mult(x,y,D):
    u,v=x
    p,q=y
    a,b=u*p+v*q*D,u*q+v*p
    assert a%2==b%2==0
    return a//2,b//2


def conj(x):
    return x[0],-x[1]


def real_sign(x,D):
    u,v=x
    if v==0:
        return (u>0)-(u<0)
    if u==0 or (u>0)==(v>0):
        return (v>0)-(v<0)
    q=u*u-v*v*D
    assert q!=0
    return ((q>0)-(q<0))*((u>0)-(u<0))


def floor_quad(x,D):
    u,v=x
    low = isqrt(v*v*D)
    if v<0:
        low=-low-1
    return (u+low)//2


def shifts(x,D):
    if x[1]==0:
        assert x[0]%2==0
        t=x[0]//2
        return ((-t,ZERO),(1-t,ONE))
    g=-floor_quad(x,D)
    y=plus(x,ci(g))
    assert real_sign(y,D)>0 and real_sign(sub(ONE,y),D)>0
    return ((g,y),)


def unit(x,D):
    return real_sign(x,D)>=0 and real_sign(sub(ONE,x),D)>=0


def product_cross(A,B,k,D):
    cross=ZERO
    for i in range(1,k):
        cross=plus(cross,mult(A[i],B[k-i],D))
    return cross


def extend(A,B,x,y,k,D,with_binary_one=False):
    cross=product_cross(A,B,k,D)
    for g,ak in shifts(x,D):
        for q,bk in shifts(y,D):
            c=plus(cross,plus(ak,bk))
            if c==ONE or (not with_binary_one and c==ZERO):
                yield g,q,ak,bk


def house(x,num,den,D):
    u,v=x
    gap=2*num-den*abs(u)
    return gap>0 and gap*gap>den*den*v*v*D


def newton_base(P,S,k,D):
    x=ZERO
    for j in range(1,k):
        x=sub(x,mult(P[j],S[k-j],D))
    return x


def coeff_base(P, G, k):
    q=ZERO
    for j in range(1,k+1):
        q=plus(q,scale(G[k-j],P[j]))
    return q


def screen_step(a):
    degree=len(a)-1
    while degree:
        leading,constant=a[0],a[degree]
        if abs(leading)<=abs(constant):
            return False
        out=[leading*a[i]-constant*a[degree-i] for i in range(degree)]
        divisor=gcd(*out)
        a=[x//divisor for x in out]
        degree-=1
    return True


def root_code(P,D):
    for p,q in ((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6)):
        value=ZERO
        for j,x in enumerate(P):
            value=plus(value,scale(p**j*q**(9-j),x))
        if real_sign(value,D)<=0 or real_sign(conj(value),D)<=0:
            return "root_positive_reject"
    H=[0]*19
    for i,(u,v) in enumerate(P):
        for j,(r,s) in enumerate(P):
            H[i+j]+=u*r-v*s*D
    assert all(h%4==0 for h in H)
    H=[h//4 for h in H]
    if not screen_step([H[j]<<j for j in range(18,-1,-1)]):
        return "root_outer_reject"
    if not screen_step([H[j]<<(18-j) for j in range(19)]):
        return "root_inner_reject"
    return "root_pass"

TOP=((0,0),(1,0),(0,1))


def forced_choices(x,D,force):
    if not force:return shifts(x,D)
    if x[1]:return ()
    assert x[0]%2==0
    return ((-x[0]//2,ZERO),)



def newton_prefix(P,D):
    active=[ci(9)]
    for k in range(1,len(P)):
        s=plus(newton_base(P,active,k,D),scale(-k,P[k]))
        limit=4*(Fraction(2)**k+Fraction(1,2)**k)+1
        if not house(s,limit.numerator,limit.denominator,D):
            return False
        active.append(s)
    return True


def prefixes(P,D,r,depth,forward):
    m,t=9-2*r,9-r
    rows=[((1,),(1,),(ONE,),(ONE,))]
    for k in range(1,depth+1):
        upcoming=[]
        for G,Q,A,B in rows:
            x,y=coeff_base(P,G,k),coeff_base(tuple(map(conj,P)),Q,k)
            for g,b,ak,bk in extend(A,B,x,y,k,D):
                if forward and k==m and (g or b):continue
                if ((forward and k==t) or (not forward and k==r)) and (ak!=ZERO or bk!=ZERO):
                    continue
                if forward and r==1 and k<=8:
                    if plus(product_cross(A,B,k,D),plus(ak,bk))!=ONE:continue
                upcoming.append((G+(g,),Q+(b,),A+(ak,),B+(bk,)))
        rows=upcoming
        if not rows:break
    return rows


def find_pair(P,D,r,depth):
    forwards=prefixes(P,D,r,depth,True)
    if not forwards:return None,"forward"
    reverses=prefixes(P[::-1],D,r,depth,False)
    if not reverses:return None,"reverse"
    for G,Q,A,B in forwards:
        for U,V,C,F in reverses:
            if all(not ((A[k]!=ZERO and F[k]!=ZERO) or (B[k]!=ZERO and C[k]!=ZERO))
                   for k in range(1,depth+1)):
                return (G,Q,A,B,U,V,C,F),"coupled"
    return None,"endpoint"


def decode(D,a,r,integer_coeffs):
    ints=iter(integer_coeffs)
    return (ONE,*(a if k==r else ci(next(ints)) for k in range(1,9)),ONE)


def run_site(r,depth):
    stats=Counter()
    witness=None
    # Only a deterministic small subbox. Every other integral coefficient and field is untested.
    for integers in ((0,1),(-1,0,1)):
        for D,trace,diff in ((5,1,1),(5,-1,1),(2,2,2),(3,2,2),(5,3,1),(13,3,1),(17,3,1),(5,0,2)):
            a=(trace,diff)
            for choices in product(integers,repeat=7):
                if len(integers)==3 and -1 not in choices:continue
                P=decode(D,a,r,choices)
                if r==4 and P[1]!=ZERO:continue
                stats["searched"]+=1
                if not newton_prefix(P,D) or not newton_prefix(P[::-1],D):
                    stats["newton_rejected"]+=1
                    continue
                pair,step=find_pair(P,D,r,depth)
                stats[step]+=1
                if pair is None:continue
                root=root_code(P,D)
                stats[root]+=1
                if root=="root_pass":
                    witness=dict(D=D,P=[[u,v] for u,v in P],r=r,depth=depth,
                                 simultaneous_prefix=pair,root_status=root)
                    return stats,witness
    return stats,witness


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--depth",type=int,default=9,choices=range(8,10))
    args=parser.parse_args()
    assert not args.out.exists()
    report={}
    for r in range(1,5):
        counts,witness=run_site(r,args.depth)
        report[r]=dict(counts=counts,witness=witness)
    args.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps({r:dict(counts=x["counts"],witness=bool(x["witness"])) for r,x in report.items()}))


if __name__=="__main__":
    main()
