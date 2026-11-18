#!/usr/bin/env python3
"""Exact eight-degree S={2,3,5} sufficient-prefix contradictions or finite entry."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path

ZERO, ONE = (0,0), (2,0)
BOUND = {2:(17,1), 3:(65,2), 4:(257,4), 5:(1025,8),
         6:(4097,16), 7:(16385,32), 8:(65537,64)}
BOX = {3:92, 4:120, 5:92, 6:41, 7:9}


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


def range_of_shift(base, k,D):
    num,den=BOUND[k]
    u,v=base
    R=(2*num-isqrt((den*v)**2*D)-1)//den
    if R<0:
        return 1,0
    return max(-BOX[k],-((R-u)//(2*k))),min(BOX[k],(u+R)//(2*k))


def newton_base(P,S,k,D):
    x=ZERO
    for j in range(1,k):
        x=sub(x,mult(P[j],S[k-j],D))
    return x


def min_d_with_bound(x,y,upper,lo,hi,D):
    """First d such that x+d*y <= upper, for y negative."""
    left,right=lo-1,hi+1
    while right-left>1:
        mid=(left+right)//2
        if real_sign(sub(plus(x,scale(mid,y)),ci(upper)),D)<=0:
            right=mid
        else:
            left=mid
    return right


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
            value=plus(value,scale(p**j*q**(8-j),x))
        if real_sign(value,D)<=0 or real_sign(conj(value),D)<=0:
            return "root_positive_reject"
    H=[0]*17
    for i,(u,v) in enumerate(P):
        for j,(r,s) in enumerate(P):
            H[i+j]+=u*r-v*s*D
    assert all(h%4==0 for h in H)
    H=[h//4 for h in H]
    if not screen_step([H[j]<<j for j in range(16,-1,-1)]):
        return "root_outer_reject"
    if not screen_step([H[j]<<(16-j) for j in range(17)]):
        return "root_inner_reject"
    return "root_pass"

TOP=((0,0),(1,0),(0,1))


def allowed_a(D,L):
    eps=int(D%4==1)
    tau=1 if eps else 2
    for r in range(1,isqrt(288//(tau*tau*D))+1):
        floor=isqrt((tau*r)**2*D)
        lo=max(-82+floor,L*L-16+floor)
        hi=min(82-floor,L*L+16-floor)
        for T in range(lo+((eps*r-lo)%2),hi+1,2):
            yield (T,tau*r),r


def top_pair(a4,b4,x,y,f):
    return (a4==ZERO or y==f) and (b4==ZERO or x==f)


def coefficients_at_top(D,L,r,u,a,a4,b4):
    eps=int(D%4==1)
    tau=1 if eps else 2
    for t in range(-r,r+1):
        if not t or r%t or t==r*u:continue
        if (40*tau*(r*u-t))**2*D>=2050**2:continue
        e=L*r//t
        if not -41<=e<=41:continue
        d0=(eps*t,tau*t)
        for v in sorted({r//t-floor_quad(d0,D),r//t-floor_quad(conj(d0),D)}):
            d=plus(d0,ci(v))
            if not house(d,185,2,D):continue
            for x,y in TOP:
                if top_pair(a4,b4,x,y,0):
                    yield d,t,e,0,x,y,False
    if not u:return
    t=r*u
    if (tau*t)**2*D>=185**2:
        return
    d0=scale(u,a)
    rho,sigma=-floor_quad(d0,D),-floor_quad(conj(d0),D)
    for x,y in TOP:
        for f in (0,1):
            if not top_pair(a4,b4,x,y,f):continue
            U1,V1=x-f,y-f
            p,q=r+L*r*U1,r+L*r*V1
            if p%t and q%t:continue
            for e in range(-41,42):
                U2,V2=-e-f*U1,-e-f*V1
                options=set()
                if p%t==0:
                    options.add(p//t+rho-e*U1-f*U2)
                if q%t==0:
                    options.add(q//t+sigma-e*V1-f*V2)
                for v in sorted(options):
                    d=plus(d0,ci(v))
                    if house(d,185,2,D):
                        yield d,t,e,f,x,y,True


def unit_choices(x,D,force):
    if not force:return shifts(x,D)
    if x[1]:return ()
    assert x[0]%2==0
    return ((-x[0]//2,ZERO),)


def reverse_coupled(P,D,U,V,A,B,low_a,low_b):
    R=P[::-1]
    S=list(map(conj,R))
    states=[(U,V,A,B)]
    for k in range(4,9):
        nxt=[]
        for U,V,A,B in states:
            cx=product_cross(A,B,k,D)
            for iu,ak in unit_choices(coeff_base(R,U,k),D,low_b[k]!=ZERO):
                for iv,bk in unit_choices(coeff_base(S,V,k),D,low_a[k]!=ZERO):
                    if plus(cx,plus(ak,bk)) in (ZERO,ONE):
                        nxt.append((U+(iu,),V+(iv,),A+(ak,),B+(bk,)))
        if not nxt:return False
        states=nxt
    return True


def one_domain(D):
    counts=Counter()
    cases={}
    for L in (-4,-3,-2,-1,1,2,3,4):
        for a,r in allowed_a(D,L):
            counts["a_cores"]+=1
            ab=conj(a)
            S2=[ci(8),ci(-L),sub(ci(L*L),scale(2,a))]
            rho,sigma=-floor_quad(a,D),-floor_quad(ab,D)
            A2,B2=plus(a,ci(rho)),plus(ab,ci(sigma))
            assert plus(A2,B2)==ONE and (rho or sigma)
            A3,B3=(ONE,ZERO,A2,ZERO),(ONE,ZERO,B2,ZERO)
            G2,Q2=(1,-L,L*L+rho),(1,-L,L*L+sigma)
            for u in range((-L**3-32+2)//3,(-L**3+32)//3+1):
                b=plus(scale(L,a),ci(u))
                if not house(b,185,2,D):continue
                counts["newton3_b"]+=1
                bb=conj(b)
                P3,Pb3=[ONE,ci(L),a,b],[ONE,ci(L),ab,bb]
                G3,Q3=G2+(-u-L*G2[2],),Q2+(-u-L*Q2[2],)
                S3=S2+[ci(-L**3-3*u)]
                ns4=newton_base(P3,S3,4,D)
                c_min,c_max=range_of_shift(ns4,4,D)
                for rho4,sigma4,a4,b4 in extend(A3,B3,coeff_base(P3+[ZERO],G3,4),
                                                coeff_base(Pb3+[ZERO],Q3,4),4,D,True):
                    counts["row4"]+=1
                    assert (a4!=ZERO)==bool(rho) and (b4!=ZERO)==bool(sigma)
                    A4,B4=A3+(a4,),B3+(b4,)
                    for d,t,e,f,x,y,exception in coefficients_at_top(D,L,r,u,a,a4,b4):
                        counts["high_options"]+=1
                        db=conj(d)
                        U1,V1=x-f,y-f
                        U2,V2=-e-f*U1,-e-f*V1
                        n,m=-floor_quad(d,D),-floor_quad(db,D)
                        U3,V3=n-e*U1-f*U2,m-e*V1-f*V2
                        force_a=r+L*r*V1+t*V3!=0
                        force_b=r+L*r*U1+t*U3!=0
                        assert not (force_a and force_b)
                        da=[]
                        if force_a:
                            if L*t%r:continue
                            da.append(rho4+L*G3[3]-L*t//r)
                        if force_b:
                            if L*t%r:continue
                            da.append(sigma4+L*Q3[3]-L*t//r)
                        if da and any(z!=da[0] for z in da):continue
                        allowed_c=(da[0],) if da else range(c_min,c_max+1)
                        AA=(ONE,ci(x),ZERO,plus(d,ci(n)))
                        BB=(ONE,ci(y),ZERO,plus(db,ci(m)))
                        U,V=(1,U1,U2,U3),(1,V1,V2,V3)
                        assert plus(AA[3],BB[3])==ONE
                        for c in allowed_c:
                            if not c_min<=c<=c_max:continue
                            P=P3+[ci(c),d,ci(e),ci(f),ONE]
                            Pb=Pb3+[ci(c),db,ci(e),ci(f),ONE]
                            S=S3+[plus(ns4,ci(-4*c))]
                            for k in range(5,9):
                                value=plus(newton_base(P,S,k,D),scale(-k,P[k]))
                                if not house(value,*BOUND[k],D):break
                                S.append(value)
                            if len(S)!=9:continue
                            counts["newton8"]+=1
                            G4,Q4=G3+(rho4-c,),Q3+(sigma4-c,)
                            for g5,q5,a5,b5 in extend(A4,B4,coeff_base(P,G4,5),coeff_base(Pb,Q4,5),5,D):
                                G5,Q5=G4+(g5,),Q4+(q5,)
                                A5,B5=A4+(a5,),B4+(b5,)
                                for g6,q6,a6,b6 in extend(A5,B5,coeff_base(P,G5,6),
                                                                coeff_base(Pb,Q5,6),6,D,True):
                                    if (force_a and a6!=ZERO) or (force_b and b6!=ZERO):continue
                                    counts["row6"]+=1
                                    G6,Q6=G5+(g6,),Q5+(q6,)
                                    A6,B6=A5+(a6,),B5+(b6,)
                                    for g7,q7,a7,b7 in extend(A6,B6,coeff_base(P,G6,7),coeff_base(Pb,Q6,7),7,D):
                                        G7,Q7=G6+(g7,),Q6+(q7,)
                                        A7,B7=A6+(a7,),B6+(b7,)
                                        for _,_,a8,b8 in extend(A7,B7,coeff_base(P,G7,8),coeff_base(Pb,Q7,8),8,D):
                                            counts["row8"]+=1
                                            if not reverse_coupled(P,D,U,V,AA,BB,A7+(a8,),B7+(b8,)):
                                                continue
                                            counts["coupled_row8"]+=1
                                            row=(D,L,a[0],r,u,c,d[0],t,e,f)
                                            if row not in cases:
                                                cases[row]=root_code(P,D)
                                                counts["mothers_exception" if exception else "mothers_generic"]+=1
    counts["mothers"]=len(cases)
    counts.update(cases.values())
    return counts,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,289) if all(D%(p*p) for p in range(2,isqrt(D)+1))
             and any(any(allowed_a(D,L)) for L in (1,2,3,4))]
    totals=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for ct,local in pool.map(one_domain,domains,chunksize=1):
            totals.update(ct)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
    payload=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in sorted(rows.items()))
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    result=dict(status="complete",domains=len(domains),counts=dict(totals),
                rows_sha256=hashlib.sha256(payload).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
