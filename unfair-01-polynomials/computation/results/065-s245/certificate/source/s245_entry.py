#!/usr/bin/env python3
"""Exact eight-degree S={2,4,5} sufficient-prefix contradictions or finite entry."""

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


def a_cores(D):
    eps=int(D%4==1)
    tau=1 if eps else 2
    for r in range(1,isqrt(288//(tau*tau*D))+1):
        R=16-isqrt((tau*r)**2*D)
        for T in range(-R+((eps*r+R)%2),R+1,2):
            yield (T,tau*r),r


def interval_minus(base,k,num,den,D):
    u,v=base
    R=(2*num-isqrt((den*v)**2*D)-1)//den
    if R<0:return 1,0
    scale_=2*k
    return -((R-u)//scale_),(u+R)//scale_


def top_options(a,r,b,t,j,a4,b4,D):
    eps=int(D%4==1)
    tau=1 if eps else 2
    exception=(t==r*b)
    if not exception:
        f=j
        if not -9<=f<=9:return
        if (40*tau*(r*b-t))**2*D>=2050**2:return
        base0=(eps*t,tau*t)
        v_options={r//t-floor_quad(base0,D),r//t-floor_quad(conj(base0),D)}
        for v in sorted(v_options):
            d=plus(base0,ci(v))
            if not house(d,185,2,D):continue
            for x,y in TOP:
                if (a4!=ZERO and y) or (b4!=ZERO and x):continue
                yield d,0,f,x,y,False
        return
    assert abs(b)==1
    d0=scale(b,a)
    qA,qB=-floor_quad(d0,D),-floor_quad(conj(d0),D)
    for x,y in TOP:
        for f in (j,j+1):
            if not -9<=f<=9:continue
            if (a4!=ZERO and y-f+j) or (b4!=ZERO and x-f+j):continue
            U1,V1=x-f,y-f
            for e in range(-41,42):
                U2,V2=-e-f*U1,-e-f*V1
                v_options={r//t+qA-e*U1+(j-f)*U2,r//t+qB-e*V1+(j-f)*V2}
                for v in sorted(v_options):
                    if not -25<=v<=25:continue
                    d=plus(d0,ci(v))
                    if house(d,185,2,D):
                        yield d,e,f,x,y,True


def forced_choices(x,D,force):
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
            cross=product_cross(A,B,k,D)
            for ug,ak in forced_choices(coeff_base(R,U,k),D,low_b[k]!=ZERO):
                for vg,bk in forced_choices(coeff_base(S,V,k),D,low_a[k]!=ZERO):
                    if plus(cross,plus(ak,bk)) in (ZERO,ONE):
                        nxt.append((U+(ug,),V+(vg,),A+(ak,),B+(bk,)))
        if not nxt:return False
        states=nxt
    return True


def one_domain(D):
    totals=Counter()
    cases={}
    eps=int(D%4==1)
    tau=1 if eps else 2
    for a,r in a_cores(D):
        totals["a_cores"]+=1
        ab=conj(a)
        g2,q2=-floor_quad(a,D),-floor_quad(ab,D)
        A2,B2=plus(a,ci(g2)),plus(ab,ci(q2))
        assert plus(A2,B2)==ONE
        A3,B3=(ONE,ZERO,A2,ZERO),(ONE,ZERO,B2,ZERO)
        bs4=scale(2,mult(a,a,D))
        for t in range(-r,r+1):
            if not t or r%t:continue
            for j in range(-10,10):
                if not j:continue
                s=t*j
                c0=(eps*s,tau*s)
                cb0=conj(c0)
                ns4=sub(bs4,scale(4,c0))
                l1,h1=interval_minus(ns4,4,257,4,D)
                l2,h2=interval_minus(neg(c0),1,1921,16,D)
                clo,chi=max(l1,l2),min(h1,h2)
                if clo>chi:continue
                totals["core_a_t_j"]+=1
                for rho4,sigma4,a4,b4 in extend(A3,B3,plus(c0,scale(g2,a)),
                                                 plus(cb0,scale(q2,ab)),4,D,True):
                    totals["row4_normalized"]+=1
                    assert (a4!=ZERO)==bool(s+r*g2) and (b4!=ZERO)==bool(s+r*q2)
                    A4,B4=A3+(a4,),B3+(b4,)
                    for b in range(-10,11):
                        P3,Pb3=[ONE,ZERO,a,ci(b)],[ONE,ZERO,ab,ci(b)]
                        S3=[ci(8),ZERO,scale(-2,a),ci(-3*b)]
                        G3,Q3=(1,0,g2,-b),(1,0,q2,-b)
                        for d,e,f,x,y,exception in top_options(a,r,b,t,j,a4,b4,D):
                            db=conj(d)
                            ns5=scale(5,sub(scale(b,a),d))
                            if not house(ns5,*BOUND[5],D):continue
                            totals["high_options"]+=1
                            U1,V1=x-f,y-f
                            U2,V2=-e-f*U1,-e-f*V1
                            rho_d,sigma_d=-floor_quad(d,D),-floor_quad(db,D)
                            U3,V3=rho_d-e*U1-f*U2,sigma_d-e*V1-f*V2
                            highA=r+s*U2+t*U3
                            highB=r+s*V2+t*V3
                            assert not (highA and highB)
                            forced=[]
                            if highB:
                                if s*g2%r:continue
                                forced.append(rho4+s*g2//r)
                            if highA:
                                if s*q2%r:continue
                                forced.append(sigma4+s*q2//r)
                            possible_c=(forced[0],) if forced else range(clo,chi+1)
                            X,Y=(ONE,ci(x),ZERO,plus(d,ci(rho_d))),(ONE,ci(y),ZERO,plus(db,ci(sigma_d)))
                            assert plus(X[3],Y[3])==ONE
                            for v_c in possible_c:
                                if not clo<=v_c<=chi:continue
                                c,cb=plus(c0,ci(v_c)),plus(cb0,ci(v_c))
                                P=P3+[c,d,ci(e),ci(f),ONE]
                                Pb=Pb3+[cb,db,ci(e),ci(f),ONE]
                                S=S3+[plus(ns4,ci(-4*v_c)),ns5]
                                for k in range(6,9):
                                    value=plus(newton_base(P,S,k,D),scale(-k,P[k]))
                                    if not house(value,*BOUND[k],D):break
                                    S.append(value)
                                if len(S)!=9:continue
                                totals["newton8"]+=1
                                G4,Q4=G3+(rho4-v_c,),Q3+(sigma4-v_c,)
                                for g5,q5,a5,b5 in extend(A4,B4,coeff_base(P,G4,5),coeff_base(Pb,Q4,5),5,D):
                                    G5,Q5=G4+(g5,),Q4+(q5,)
                                    A5,B5=A4+(a5,),B4+(b5,)
                                    for g6,q6,a6,b6 in extend(A5,B5,coeff_base(P,G5,6),
                                                             coeff_base(Pb,Q5,6),6,D,True):
                                        if (highB and a6!=ZERO) or (highA and b6!=ZERO):continue
                                        totals["row6"]+=1
                                        G6,Q6=G5+(g6,),Q5+(q6,)
                                        A6,B6=A5+(a6,),B5+(b6,)
                                        for g7,q7,a7,b7 in extend(A6,B6,coeff_base(P,G6,7),
                                                                 coeff_base(Pb,Q6,7),7,D):
                                            A7,B7=A6+(a7,),B6+(b7,)
                                            G7,Q7=G6+(g7,),Q6+(q7,)
                                            for _,_,a8,b8 in extend(A7,B7,coeff_base(P,G7,8),coeff_base(Pb,Q7,8),8,D):
                                                totals["row8"]+=1
                                                if not reverse_coupled(P,D,(1,U1,U2,U3),(1,V1,V2,V3),X,Y,
                                                                       A7+(a8,),B7+(b8,)):
                                                    continue
                                                totals["coupled_row8"]+=1
                                                row=(D,a[0],r,b,c[0],s,d[0],t,e,f)
                                                if row not in cases:
                                                    cases[row]=root_code(P,D)
                                                    totals["mothers_exception" if exception else "mothers_generic"]+=1
    totals["mothers"]=len(cases)
    totals.update(cases.values())
    return totals,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,289) if all(D%(p*p) for p in range(2,isqrt(D)+1))
             and any(a_cores(D))]
    counts=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for ct,local in pool.map(one_domain,domains,chunksize=1):
            counts.update(ct)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
    payload=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in sorted(rows.items()))
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    result=dict(status="complete",domains=len(domains),counts=dict(counts),
                rows_sha256=hashlib.sha256(payload).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
