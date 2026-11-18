#!/usr/bin/env python3
"""Exact eight-degree S={3,4} sufficient-prefix contradictions or finite entry."""

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


def allowed_a(D,L,M):
    t=int(D%4==1)
    tau=1 if t else 2
    vmax=isqrt((65**2-1)//(9*tau*tau*D))
    base=-2*L**3+6*L*M
    for r in range(1,vmax+1):
        radical=isqrt((3*tau*r)**2*D)
        margin=64-radical
        hlo=(base-margin+2)//3
        hhi=(base+margin)//3
        amax=184-isqrt((tau*r)**2*D)
        hlo=max(hlo,-amax)
        hhi=min(hhi,amax)
        for T in range(hlo+((t*r-hlo)%2),hhi+1,2):
            yield (T,tau*r),r


def e_d_choices(L,A2,B2,A5,B5,A6,B6):
    for x,y in TOP:
        for p,q in TOP:
            if A2*q or B2*p:continue
            for e in range(-9,10):
                u1,v1=x-e,y-e
                if (B5!=ZERO and 1+L*u1) or (A5!=ZERO and 1+L*v1):
                    continue
                ds=[]
                if B6!=ZERO:
                    dA=L*p+(1-L*e)*u1
                    if dA%L:continue
                    ds.append(dA//L)
                if A6!=ZERO:
                    dB=L*q+(1-L*e)*v1
                    if dB%L:continue
                    ds.append(dB//L)
                assert ds
                d=ds[0]
                if not -41<=d<=41 or any(i!=d for i in ds):
                    continue
                u2,v2=p-d-e*u1,q-d-e*v1
                yield x,y,p,q,e,d,u1,v1,u2,v2


def c_choices(L,M,u,g2,q2,d,e,u1,v1,u2,v2,clo,chi):
    needed=[]
    if u+M*g2:
        if v2%L:return ()
        needed.append(v2//L-d*v1-e*v2)
    if u+M*q2:
        if u2%L:return ()
        needed.append(u2//L-d*u1-e*u2)
    if not needed:
        return range(clo,chi+1)
    c=needed[0]
    return (c,) if clo<=c<=chi and all(y==c for y in needed) else ()


def force(base,D,needed):
    if not needed:
        return shifts(base,D)
    if base[1]:
        return ()
    return ((-base[0]//2,ZERO),)


def reverse_coupled(P,D,G,Q,A,B,low_a,low_b):
    R=P[::-1]
    S=list(map(conj,R))
    states=[(G,Q,A,B)]
    for k in range(4,9):
        nxt=[]
        for U,V,C,E in states:
            cx=product_cross(C,E,k,D)
            for ui,ck in force(coeff_base(R,U,k),D,low_b[k]!=ZERO):
                for vi,ek in force(coeff_base(S,V,k),D,low_a[k]!=ZERO):
                    if plus(cx,plus(ck,ek)) in (ZERO,ONE):
                        nxt.append((U+(ui,),V+(vi,),C+(ck,),E+(ek,)))
        if not nxt:
            return False
        states=nxt
    return True


def one_domain(D):
    counts=Counter()
    cases={}
    for L in (-4,-3,-2,-1,1,2,3,4):
        for M in (range(-7,9) if abs(L)==1 else (0,)):
            for a,r in allowed_a(D,L,M):
                counts["a_cores"]+=1
                ab=conj(a)
                s2=ci(L*L-2*M)
                s3=plus(ci(-L**3+3*L*M),scale(-3,a))
                assert house(a,185,2,D) and house(s3,*BOUND[3],D)
                s4_const=L**4-4*L*L*M+2*M*M
                u_lo=(s4_const-64+3)//4
                u_hi=(s4_const+64)//4
                rho,sigma=-floor_quad(a,D),-floor_quad(ab,D)
                A3,B3=plus(a,ci(rho)),plus(ab,ci(sigma))
                assert plus(A3,B3)==ONE
                S3=[ci(8),ci(-L),s2,s3]
                for a2,b2 in TOP if abs(L)==1 else ((0,0),):
                    g2=L*L-M+a2
                    q2=L*L-M+b2
                    g3=rho+L*M-L*g2
                    q3=sigma+L*M-L*q2
                    G3,Q3=(1,-L,g2,g3),(1,-L,q2,q3)
                    P3=[ONE,ci(L),ci(M),a]
                    Pb3=list(map(conj,P3))
                    A4,B4=(ONE,ZERO,ci(a2),A3,ZERO),(ONE,ZERO,ci(b2),B3,ZERO)
                    P4_temp,Pb4_temp=P3+[scale(L,a)],Pb3+[scale(L,ab)]
                    g4_const=-M*g2-L*g3
                    q4_const=-M*q2-L*q3
                    G4t,Q4t=G3+(g4_const,),Q3+(q4_const,)
                    P5t,Pb5t=P4_temp+[ZERO],Pb4_temp+[ZERO]
                    for g5t,q5t,a5,b5 in extend(A4,B4,coeff_base(P5t,G4t,5),coeff_base(Pb5t,Q4t,5),5,D):
                        counts["row5_core"]+=1
                        A5,B5=A4+(a5,),B4+(b5,)
                        G5t,Q5t=G4t+(g5t,),Q4t+(q5t,)
                        for _,_,a6,b6 in extend(A5,B5,coeff_base(P5t+[ZERO],G5t,6),
                                                 coeff_base(Pb5t+[ZERO],Q5t,6),6,D,True):
                            counts["row6_core"]+=1
                            assert (a6==ZERO)==(rho+M*L==0)
                            assert (b6==ZERO)==(sigma+M*L==0)
                            A6,B6=A5+(a6,),B5+(b6,)
                            for x,y,p,q,e,d,u1,v1,u2,v2 in e_d_choices(L,a2,b2,a5,b5,a6,b6):
                                counts["top_options"]+=1
                                for u in range(u_lo,u_hi+1):
                                    b=plus(scale(L,a),ci(u))
                                    bb=conj(b)
                                    if not house(b,1921,16,D):
                                        continue
                                    P4,Pb4=P3+[b],Pb3+[bb]
                                    S4=S3+[ci(s4_const-4*u)]
                                    ns5=newton_base(P4,S4,5,D)
                                    clo,chi=range_of_shift(ns5,5,D)
                                    for c in c_choices(L,M,u,g2,q2,d,e,u1,v1,u2,v2,clo,chi):
                                        counts["newton5_top_c"]+=1
                                        P=P4+[ci(c),ci(d),ci(e),ONE]
                                        Pb=Pb4+[ci(c),ci(d),ci(e),ONE]
                                        S=S4+[plus(ns5,ci(-5*c))]
                                        for k,co in ((6,d),(7,e),(8,1)):
                                            v=plus(newton_base(P,S,k,D),ci(-k*co))
                                            if not house(v,*BOUND[k],D):break
                                            S.append(v)
                                        if len(S)!=9:continue
                                        counts["newton8_top_c"]+=1
                                        G4,Q4=G3+(g4_const-u,),Q3+(q4_const-u,)
                                        ga5,gb5=sub(a5,coeff_base(P,G4,5)),sub(b5,coeff_base(Pb,Q4,5))
                                        assert ga5[1]==gb5[1]==0 and ga5[0]%2==gb5[0]%2==0
                                        G5,Q5=G4+(ga5[0]//2,),Q4+(gb5[0]//2,)
                                        ga6,gb6=sub(a6,coeff_base(P,G5,6)),sub(b6,coeff_base(Pb,Q5,6))
                                        assert ga6[1]==gb6[1]==0 and ga6[0]%2==gb6[0]%2==0
                                        G6,Q6=G5+(ga6[0]//2,),Q5+(gb6[0]//2,)
                                        for g7,q7,a7,b7 in extend(A6,B6,coeff_base(P,G6,7),coeff_base(Pb,Q6,7),7,D):
                                            A7,B7=A6+(a7,),B6+(b7,)
                                            G7,Q7=G6+(g7,),Q6+(q7,)
                                            counts["row7"]+=1
                                            for _,_,a8,b8 in extend(A7,B7,coeff_base(P,G7,8),coeff_base(Pb,Q7,8),8,D):
                                                counts["row8"]+=1
                                                U3,V3=-c-d*u1-e*u2,-c-d*v1-e*v2
                                                if not reverse_coupled(P,D,(1,u1,u2,U3),(1,v1,v2,V3),
                                                    (ONE,ci(x),ci(p),ZERO),(ONE,ci(y),ci(q),ZERO),A7+(a8,),B7+(b8,)):
                                                    continue
                                                counts["coupled_row8"]+=1
                                                key=(D,L,M,a[0],r,u,c,d,e)
                                                if key not in cases:
                                                    cases[key]=root_code(P,D)
    counts["mothers"]=len(cases)
    counts.update(cases.values())
    return counts,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,470) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    totals=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for counts,local in pool.map(one_domain,domains,chunksize=1):
            totals.update(counts)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
    raw=b"".join(json.dumps(k,separators=(",",":")).encode()+b"\n" for k in sorted(rows.items()))
    a.out.with_suffix(".rows.jsonl").write_bytes(raw)
    out=dict(status="complete",domains_scanned=len(domains),counts=dict(totals),rows_sha256=hashlib.sha256(raw).hexdigest())
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__=="__main__":
    main()
