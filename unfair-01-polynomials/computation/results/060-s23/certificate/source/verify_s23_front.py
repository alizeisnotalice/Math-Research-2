#!/usr/bin/env python3
"""Independent integral-basis S23 replay using complete integral-basis coordinates and direct common shifts."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path

Z=(0,0)
I=(1,0)
BOUND={2:(17,1),3:(65,2),4:(257,4),5:(1025,8),6:(4097,16),7:(16385,32),8:(65537,64)}


def add(x,y): return x[0]+y[0],x[1]+y[1]
def mulc(k,x): return k*x[0],k*x[1]
def integer(x): return x,0


def audit_domain(D):
    t=int(D%4==1)
    radicand_factor=1 if t else 2
    v=(D-1)//4 if t else D

    def mul(x,y):
        p,q=x
        r,s=y
        return p*r+q*s*v,p*s+q*r+q*s*t

    def conjugate(x):
        return x[0]+t*x[1],-x[1]

    def sign(x):
        m,n=x
        u=2*m+t*n
        w=radicand_factor*n
        if u>=0 and w>=0: return int(u+w>0)
        if u<=0 and w<=0: return -int(u+w<0)
        q=u*u-D*w*w
        assert q
        return int(q>0)*2-1 if u>0 else 1-int(q>0)*2

    def both_abs(x,num,den):
        m,n=x
        gap=2*num-den*abs(2*m+t*n)
        return gap>0 and gap*gap>D*(den*radicand_factor*n)**2

    def choices(x):
        m,n=x
        if not n:
            return ((-m,Z),(1-m,I))
        R=abs(m)+abs(n)*(D+1)+2
        lo,hi=-R,R
        while hi-lo>1:
            mid=(lo+hi)//2
            if sign((m-mid,n))>=0: lo=mid
            else: hi=mid
        return ((-lo,(m-lo,n)),)

    def engine(P,G,k):
        out=Z
        for i in range(1,k+1):
            out=add(out,mulc(G[k-i],P[i]))
        return out

    def cross(A,B,k):
        out=Z
        for i in range(1,k):
            out=add(out,mul(A[i],B[k-i]))
        return out

    def update(A,B,u,w,k,req=(Z,I)):
        h=cross(A,B,k)
        return [(g,q,a,b) for g,a in choices(u) for q,b in choices(w)
                if add(add(a,b),h) in req]

    def newton_base(P,S,k):
        base=Z
        for j in range(1,k):
            base=add(base,mulc(-1,mul(P[j],S[k-j])))
        return base

    counts=Counter()
    rows=[]
    for r in range(1,9):
        for m in range(-45,46):
            a=(m,r)
            if not both_abs(a,83,2):
                continue
            ab=conjugate(a)
            (ka,a2),=choices(a)
            (kb,b2),=choices(ab)
            assert add(a2,b2)==I
            for L in range(-9,10):
                s2=add(integer(L*L),mulc(-2,a))
                if not both_abs(s2,17,1):
                    continue
                counts["a_cores"]+=1
                k3=add(integer(-L**3),mulc(3*L,a))
                for A1,B1 in ((0,0),(1,0),(0,1)):
                    g1,q1=A1-L,B1-L
                    g2,q2=ka-L*g1,kb-L*q1
                    A,B=(I,integer(A1),a2),(I,integer(B1),b2)
                    for s in range(-84,85):
                        if not s or D*(radicand_factor*s)**2>=185**2:
                            continue
                        if D*(3*radicand_factor*(L*r-s))**2>=65**2:
                            continue
                        counts["irr_core"]+=1
                        if s*g2%r or s*q2%r:
                            continue
                        counts["row5_projection_divisibility"]+=1
                        g3,q3=-s*g2//r,-s*q2//r
                        bstart,bbarstart=(0,s),(t*s,-s)
                        left=add(bstart,add(mulc(g1,a),integer(L*g2+g3)))
                        right=add(bbarstart,add(mulc(q1,ab),integer(L*q2+q3)))
                        internal=cross(A,B,3)
                        for u,a3 in choices(left):
                            b3=add(right,integer(u))
                            if sign(b3)<0 or sign(add(I,mulc(-1,b3)))<0:
                                continue
                            if add(add(a3,b3),internal) not in (Z,I):
                                continue
                            counts["row3_unit"]+=1
                            b=(u,s)
                            bb=conjugate(b)
                            s3=add(k3,mulc(-3,b))
                            if not both_abs(b,185,2) or not both_abs(s3,65,2):
                                continue
                            counts["row3"]+=1
                            A3,B3=A+(a3,),B+(b3,)
                            P3,PB3=[I,integer(L),a,b],[I,integer(L),ab,bb]
                            G3,Q3=(1,g1,g2,g3),(1,q1,q2,q3)
                            for rho4,sigma4,a4,b4 in update(
                                    A3,B3,engine(P3+[Z],G3,4),engine(PB3+[Z],Q3,4),4,(I,)):
                                counts["row4_normalized"]+=1
                                budget=cross(A3+(a4,),B3+(b4,),5)
                                assert budget[1]==0
                                counts[f"C5_{budget[0]}"]+=1
                                if budget not in (Z,I):
                                    continue
                                counts["row5"]+=1
                                counts[f"row5_low{A1}{B1}"]+=1
                                rows.append((D,L,2*m+t*r,r,2*u+t*s,s,A1,B1,
                                             *G3[1:],*Q3[1:],rho4,sigma4))
    return counts,rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,289) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    totals=Counter()
    rows=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for counts,local in pool.map(audit_domain,domains,chunksize=1):
            totals.update(counts)
            rows.extend(local)
    rows.sort()
    payload=b"".join(json.dumps(r,separators=(",",":")).encode()+b"\n" for r in rows)
    assert payload==a.reference.read_bytes()
    expected=json.loads(a.reference.with_name("front.json").read_bytes())
    assert totals.pop("a_cores")==expected["a_cores"] and dict(totals)==expected["counts"]
    result=dict(status="verified",a_cores=expected["a_cores"],counts=dict(totals),
                rows_sha256=sha256(payload).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
