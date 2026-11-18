#!/usr/bin/env python3
"""Exact eight-degree S={1,2,3} sufficient-prefix contradictions or finite entry."""

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


TOP=((0,0),(1,0),(0,1))


def interval_minus(base,k,num,den,D):
    u,v=base
    radius=(2*num-isqrt((den*v)**2*D)-1)//den
    if radius<0:return 1,0
    scale_=2*k
    return -((radius-u)//scale_),(u+radius)//scale_


def a_cores(D):
    tau=1 if D%4==1 else 2
    t=int(D%4==1)
    for r in range(1,isqrt(399//(tau*tau*D))+1):
        R=19-isqrt((tau*r)**2*D)
        for T in range(-R+((t*r+R)%2),R+1,2):
            yield (T,r*tau),r


def b_coordinates(a,D,g):
    tau=1 if D%4==1 else 2
    eps=int(D%4==1)
    aa=mult(a,a,D)
    tbound=isqrt((83**2-1)//(tau*tau*D))//abs(g)
    for t in range(-tbound,tbound+1):
        if not t:continue
        s=t*g
        rad=isqrt((aa[1]-2*tau*s)**2*D)
        R=33-rad
        if R<0:continue
        house_b=82-isqrt((tau*s)**2*D)
        lo=max(-house_b,-((R-aa[0])//2))
        hi=min(house_b,(R+aa[0])//2)
        for T in range(lo+((eps*s-lo)%2),hi+1,2):
            yield (T,tau*s),s,t


def one_domain(D):
    counts=Counter()
    rows={}
    tau=1 if D%4==1 else 2
    eps=int(D%4==1)
    for a,r in a_cores(D):
        counts["a_cores"]+=1
        ab=conj(a)
        g1,q1=-floor_quad(a,D),-floor_quad(ab,D)
        A=(ONE,plus(a,ci(g1)))
        B=(ONE,plus(ab,ci(q1)))
        assert plus(A[1],B[1])==ONE
        aa=mult(a,a,D)
        slope=sub(ci(g1-q1),sub(a,ab))
        assert real_sign(plus(slope,ci(2)),D)<0
        for g in range(-9,10):
            if g==0:continue
            for b,s,t in b_coordinates(a,D,g):
                counts["newton2_ab"]+=1
                bb=conj(b)
                s2=sub(aa,scale(2,b))
                assert house(b,83,2,D) and house(s2,17,1,D)
                for g2,q2,a2,b2 in extend(A,B,plus(b,scale(g1,a)),plus(bb,scale(q1,ab)),2,D,True):
                    counts["row2"]+=1
                    A2,B2=A+(a2,),B+(b2,)
                    c0=(eps*t,tau*t)
                    cb0=conj(c0)
                    U3=plus(c0,plus(scale(g1,b),scale(g2,a)))
                    V3=plus(cb0,plus(scale(q1,bb),scale(q2,ab)))
                    for rho3,sigma3,a3,b3 in extend(A2,B2,U3,V3,3,D,True):
                        counts["row3_normalized"]+=1
                        R=s*(rho3-sigma3)+t*(g2-q2)
                        if R%r:continue
                        counts["row5_difference_divisibility"]+=1
                        gap0=plus(ci(-R//r),
                                  plus(sub(scale(rho3,a),scale(sigma3,ab)),
                                       plus(sub(scale(g2,b),scale(q2,bb)),
                                            sub(scale(g1,c0),scale(q1,cb0)))))
                        P2=[ONE,a,b]
                        S2=[ci(8),neg(a),s2]
                        ns3=plus(newton_base(P2,S2,3,D),scale(-3,c0))
                        v_min1,v_max1=interval_minus(ns3,3,65,2,D)
                        v_min2,v_max2=interval_minus(neg(c0),1,185,2,D)
                        lo,hi=max(v_min1,v_min2),min(v_max1,v_max2)
                        if lo>hi:continue
                        v=min_d_with_bound(gap0,slope,1,lo,hi,D)
                        if v>hi or real_sign(plus(plus(gap0,scale(v,slope)),ci(1)),D)<0:continue
                        counts["row4_difference"]+=1
                        c=plus(c0,ci(v))
                        cb=conj(c)
                        g3,q3=rho3-v,sigma3-v
                        gg,qq=s*g3+t*g2,s*q3+t*q2
                        if gg%r or qq%r:continue
                        counts["row5_full_divisibility"]+=1
                        g4,q4=-gg//r,-qq//r
                        P3,Pb3=P2+[c],[ONE,ab,bb,cb]
                        G3,Q3=(1,g1,g2,g3),(1,q1,q2,q3)
                        U4,V4=plus(coeff_base(P3+[ZERO],G3,4),ci(g4)),plus(coeff_base(Pb3+[ZERO],Q3,4),ci(q4))
                        A3,B3=A2+(a3,),B2+(b3,)
                        K=plus(product_cross(A3,B3,4,D),plus(U4,V4))
                        assert K[1]==0 and K[0]%2==0
                        N=K[0]//2
                        if (1-N)%2:continue
                        counts["row4_parity"]+=1
                        d=(1-N)//2
                        a4,b4=plus(U4,ci(d)),plus(V4,ci(d))
                        if not unit(a4,D) or not unit(b4,D):continue
                        assert plus(product_cross(A3,B3,4,D),plus(a4,b4))==ONE
                        counts["shared_d_binary"]+=1
                        if not -120<=d<=120:continue
                        S3=S2+[plus(ns3,ci(-3*v))]
                        s4=plus(newton_base(P3,S3,4,D),ci(-4*d))
                        if not house(s4,*BOUND[4],D):continue
                        counts["row4"]+=1
    counts["mothers"]=len(rows)
    counts.update(rows.values())
    return counts,rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,400) if all(D%(p*p) for p in range(2,isqrt(D)+1))
             and any(a_cores(D))]
    totals=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for cs,res in pool.map(one_domain,domains,chunksize=1):
            totals.update(cs)
            assert rows.keys().isdisjoint(res)
            rows.update(res)
    assert totals["shared_d_binary"]==0
    payload=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in sorted(rows.items()))
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    result=dict(status="complete",domains=len(domains),counts=dict(totals),
                rows_sha256=hashlib.sha256(payload).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
