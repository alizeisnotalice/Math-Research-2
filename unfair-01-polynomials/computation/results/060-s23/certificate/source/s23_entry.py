#!/usr/bin/env python3
"""Exact eight-degree S={2,3} sufficient-prefix contradictions or finite entry."""

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
 
def cores():
    radicands=[D for D in range(2,289) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    for D in radicands:
        tau=1 if D%4==1 else 2
        eps=1 if D%4==1 else 0
        for r in range(1,isqrt(288//(tau*tau*D))+1):
            gap=isqrt(D*(tau*r)**2)
            for L in range(-9,10):
                hlo=max(L*L-16+gap,-82+gap)
                hhi=min(L*L+16-gap,82-gap)
                for T in range(hlo+((eps*r-hlo)%2),hhi+1,2):
                    yield D,L,T,r


def prefixes(job):
    D,L,T,r=job
    eps=1 if D%4==1 else 0
    tau=1 if D%4==1 else 2
    a=(T,tau*r)
    ab=conj(a)
    s2=sub(ci(L*L),scale(2,a))
    assert house(s2,17,1,D) and house(a,83,2,D)
    base3=plus(ci(-L**3),scale(3*L,a))
    sa_max=isqrt((185*185-1)//(tau*tau*D))
    sa_tol=isqrt((65*65-1)//(9*tau*tau*D))
    sa_min=max(-sa_max,L*r-sa_tol)
    sa_hi=min(sa_max,L*r+sa_tol)
    ak=-floor_quad(a,D)
    bk=1-T-ak
    t,w=plus(a,ci(ak)),plus(ab,ci(bk))
    counts=Counter()
    rows=[]
    for A1,B1 in ((0,0),(1,0),(0,1)):
        g1,q1=A1-L,B1-L
        g2,q2=ak-L*g1,bk-L*q1
        A,B=(ONE,ci(A1),t),(ONE,ci(B1),w)
        cross3=product_cross(A,B,3,D)
        for s in range(sa_min,sa_hi+1):
            if not s: continue
            counts["irr_core"]+=1
            if s*g2%r or s*q2%r: continue
            counts["row5_projection_divisibility"]+=1
            g3,q3=-s*g2//r,-s*q2//r
            p,q=(eps*s,tau*s),(eps*s,-tau*s)
            x3=plus(p,plus(scale(g1,a),ci(L*g2+g3)))
            y3=plus(q,plus(scale(q1,ab),ci(L*q2+q3)))
            K=plus(cross3,plus(x3,y3))
            assert K[1]==0 and K[0]%2==0
            u=-(K[0]//4)
            a3,b3=plus(x3,ci(u)),plus(y3,ci(u))
            if not unit(a3,D) or not unit(b3,D): continue
            assert plus(K,ci(2*u)) in (ZERO,ONE)
            counts["row3_unit"]+=1
            b=plus(p,ci(u))
            bb=conj(b)
            s3=plus(base3,scale(-3,b))
            if not house(b,185,2,D) or not house(s3,*BOUND[3],D): continue
            counts["row3"]+=1
            A3,B3=A+(a3,),B+(b3,)
            G3,Q3=(1,g1,g2,g3),(1,q1,q2,q3)
            P3,Pb3=[ONE,ci(L),a,b],[ONE,ci(L),ab,bb]
            row4=list(extend(A3,B3,coeff_base(P3+[ZERO],G3,4),
                             coeff_base(Pb3+[ZERO],Q3,4),4,D,True))
            assert len(row4)<=1
            for g4base,q4base,a4,b4 in row4:
                counts["row4_normalized"]+=1
                A4,B4=A3+(a4,),B3+(b4,)
                C5=product_cross(A4,B4,5,D)
                assert C5[1]==0 and C5[0]%2==0
                counts[f"C5_{C5[0]//2}"]+=1
                if C5 not in (ZERO,ONE):continue
                counts["row5"]+=1
                counts[f"row5_low{A1}{B1}"]+=1
                record=[D,L,T,r,b[0],s,A1,B1,*G3[1:],*Q3[1:],g4base,q4base]
                rows.append(record)
    return counts,rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domain=list(cores())
    counts=Counter()
    records=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for cs,rows in pool.map(prefixes, domain, chunksize=16):
            counts.update(cs)
            records.extend(rows)
    records.sort()
    assert len(records)==len(set(map(tuple,records)))==counts["row5"]
    raw=b"".join(json.dumps(row,separators=(",",":")).encode()+b"\n" for row in records)
    a.out.with_suffix(".rows.jsonl").write_bytes(raw)
    out=dict(status="complete",a_cores=len(domain),domains=len({x[0] for x in domain}),
             counts=dict(counts),rows_sha256=hashlib.sha256(raw).hexdigest())
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__=="__main__":
    main()


