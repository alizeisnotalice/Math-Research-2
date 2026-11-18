#!/usr/bin/env python3
"""Exact eight-degree S={1,3} sufficient-prefix contradictions or finite entry."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import isqrt
from pathlib import Path

ZERO, ONE = (0,0), (2,0)
BOUND = {2:(17,1), 3:(65,2), 4:(257,4), 5:(1025,8),
         6:(4097,16), 7:(16385,32), 8:(65537,64)}
BOX = {2:41, 3:92, 4:120, 5:92, 6:41, 7:9}


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


def coeff_base(P, G, k):
    q=ZERO
    for j in range(1,k+1):
        q=plus(q,scale(G[k-j],P[j]))
    return q


def candidates(D,house_two):
    tau=1 if D%4==1 else 2
    parity=1 if D%4==1 else 0
    house_double=83 if house_two else 20
    max_coef=isqrt((house_double**2-1)//(tau*tau*D))
    for n in range(-max_coef if house_two else 1,max_coef+1):
        if n==0:
            continue
        bound=house_double-isqrt((tau*n)**2*D)-1
        for T in range(-bound+((-bound-parity*n)%2),bound+1,2):
            assert (T-parity*n)%2==0
            yield (T,tau*n),n

def b_intervals(D,base):
    tau=1 if D%4==1 else 2
    eps=1 if D%4==1 else 0
    max_b=isqrt((185**2-1)//(tau*tau*D))
    u,v=base
    for s in range(-max_b,max_b+1):
        if s==0:
            continue
        b_box=184-isqrt((tau*s)**2*D)
        R=(130-isqrt((2*(v-3*tau*s))**2*D)-1)//2
        if R<0:
            continue
        lo=max(-b_box,(u-R+2)//3)
        hi=min(b_box,(u+R)//3)
        start=lo+((eps*s-lo)%2)
        for T in range(start,hi+1,2):
            yield (T,tau*s),s


def one_domain(D):
    counts=Counter()
    rows={}
    for a,r in candidates(D,False):
        counts["a_cores"]+=1
        ab=conj(a)
        g1=-floor_quad(a,D)
        q1=1-a[0]-g1
        A,B=(ONE,plus(a,ci(g1))),(ONE,plus(ab,ci(q1)))
        assert plus(A[1],B[1])==ONE
        row2=list(extend(A,B,scale(g1,a),scale(q1,ab),2,D,True))
        assert len(row2)<=1
        if not row2:
            continue
        counts["a_row2"]+=1
        rho2,sigma2,A2,B2=row2[0]
        A,B=A+(A2,),B+(B2,)
        P1=[ONE,a]
        S1=[ci(8),neg(a)]
        base2=newton_base(P1,S1,2,D)
        clo,chi=range_of_shift(base2,2,D)
        for c in range(clo,chi+1):
            counts["newton2_c"]+=1
            g2,q2=rho2-c,sigma2-c
            P2=P1+[ci(c)]
            S2=S1+[plus(base2,ci(-2*c))]
            base3=newton_base(P2,S2,3,D)
            for b,s in b_intervals(D,base3):
                s3=plus(base3,scale(-3,b))
                assert house(b,185,2,D) and house(s3,*BOUND[3],D)
                counts["newton3_core"]+=1
                if s*g2%r or s*q2%r:
                    continue
                counts["row5_projection_divisibility"]+=1
                P3=P2+[b]
                Pb3=list(map(conj,P3))
                G2,Q2=(1,g1,g2),(1,q1,q2)
                row3=list(extend(A,B,coeff_base(P3,G2,3),coeff_base(Pb3,Q2,3),3,D,True))
                assert len(row3)<=1
                if not row3:
                    continue
                counts["row3"]+=1
                g3,q3,A3,B3=row3[0]
                Athree,Bthree=A+(A3,),B+(B3,)
                G3,Q3=G2+(g3,),Q2+(q3,)
                row4=list(extend(Athree,Bthree,coeff_base(P3+[ZERO],G3,4),
                                 coeff_base(Pb3+[ZERO],Q3,4),4,D,True))
                assert len(row4)<=1
                if not row4:
                    continue
                counts["row4_normalized"]+=1
                rho4,sigma4,A4,B4=row4[0]
                g4,q4=-s*g2//r,-s*q2//r
                d=rho4-g4
                if d==sigma4-q4:
                    counts["shared_d"]+=1
                    if -120<=d<=120:
                        counts["shared_d_box"]+=1
                        S3=S2+[s3]
                        s4=plus(newton_base(P3,S3,4,D),ci(-4*d))
                        if house(s4,*BOUND[4],D):
                            counts["row4"]+=1
    counts["mothers"]=0
    return counts,rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    args=ap.parse_args()
    assert not args.out.exists()
    domains=[D for D in range(2,400) if all(D%(p*p)!=0 for p in range(2,isqrt(D)+1))
             and any(candidates(D,False))]
    all_rows={}
    counts=Counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for local,rows in pool.map(one_domain,domains,chunksize=1):
            counts.update(local)
            assert all_rows.keys().isdisjoint(rows)
            all_rows.update(rows)
    assert counts["shared_d"]==0
    cases=sorted((row,code) for row,code in all_rows.items())
    payload=b"".join(json.dumps(rec,separators=(",",":")).encode()+b"\n" for rec in cases)
    args.out.with_suffix(".rows.jsonl").write_bytes(payload)
    summary=dict(status="complete",domains=len(domains),counts=dict(counts),
                 cases_sha256=hashlib.sha256(payload).hexdigest(),
                 initial_domain_range=[min(domains),max(domains)])
    args.out.write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print(json.dumps(summary,sort_keys=True))


if __name__=="__main__":
    main()

