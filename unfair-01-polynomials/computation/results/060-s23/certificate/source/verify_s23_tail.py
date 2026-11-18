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


def basis(D):
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

    def box(base,k,M):
        m,n=base
        num,den={4:(257,4),5:(1025,8),6:(4097,16)}[k]
        R=(2*num-isqrt(D*(den*radicand_factor*n)**2)-1)//den
        if R<0:
            return range(0)
        u=2*m+t*n
        return range(max(-M,(u-R+2*k-1)//(2*k)),min(M,(u+R)//(2*k))+1)

    def annulus(P):
        candidates=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))
        for i,(p,q) in enumerate(candidates):
            evaluation=Z
            for j,x in enumerate(P):
                evaluation=add(evaluation,mulc(p**j*q**(8-j),x))
            if min(sign(evaluation),sign(conjugate(evaluation)))<=0:
                return 10+i
        H=[Z]*17
        Pb=list(map(conjugate,P))
        for i,u in enumerate(P):
            for j,v in enumerate(Pb):
                H[i+j]=add(H[i+j],mul(u,v))
        assert all(x[1]==0 for x in H)
        h=[x[0] for x in H]

        def schur(a):
            for n in range(16,0,-1):
                first,last=a[0],a[-1]
                if abs(first)<=abs(last):
                    return n
                v=[first*a[i]-last*a[n-i] for i in range(n)]
                factor=gcd(*v)
                a=[x//factor for x in v]
            return 0
        a=schur([h[j]<<j for j in range(16,-1,-1)])
        if a:
            return 30+a
        a=schur([h[j]<<(16-j) for j in range(17)])
        return 50+a if a else 255

    return locals()


TOP=((0,0),(1,0),(0,1))


def solve(row):
    D,L,Ta,r,Tb,s,A1,B1,g1,g2,g3,q1,q2,q3,rho4,sigma4=row
    env=basis(D)
    t=env["t"]
    mul,sign,conj=env["mul"],env["sign"],env["conjugate"]
    update,engine,base=env["update"],env["engine"],env["newton_base"]
    cross,box,house=env["cross"],env["box"],env["both_abs"]
    a,b=((Ta-t*r)//2,r),((Tb-t*s)//2,s)
    P3=[I,integer(L),a,b]
    Pb3=list(map(conj,P3))
    G3,Q3=(1,g1,g2,g3),(1,q1,q2,q3)
    A=[I]
    B=[I]
    for k in range(1,4):
        A.append(add(engine(P3,G3,k),integer(G3[k])))
        B.append(add(engine(Pb3,Q3,k),integer(Q3[k])))
    A.extend((add(engine(P3+[Z],G3,4),integer(rho4)),Z))
    B.extend((add(engine(Pb3+[Z],Q3,4),integer(sigma4)),Z))
    A,B=tuple(A),tuple(B)
    S3=[integer(8)]
    for k in range(1,4):
        S3.append(add(base(P3,S3,k),mulc(-k,P3[k])))
    ns4=base(P3,S3,4)
    found={}
    counters=Counter()
    for c in box(ns4,4,120):
        G4,Q4=G3+(rho4-c,),Q3+(sigma4-c,)
        options=[]
        for ea,eb in TOP:
            if A1*eb or B1*ea:
                continue
            for f in range(-9,10):
                if r+s*(eb-f) and r*G4[4]+s*G4[3]:continue
                if r+s*(ea-f) and r*Q4[4]+s*Q4[3]:continue
                options.append((f,ea,eb))
        if not options:
            continue
        counters["c_groups"]+=1
        P4,Pb4=P3+[integer(c)],Pb3+[integer(c)]
        S4=S3+[add(ns4,integer(-4*c))]
        ns5=base(P4,S4,5)
        for d in box(ns5,5,92):
            counters["newton5_cd"]+=1
            P5,Pb5=P4+[integer(d)],Pb4+[integer(d)]
            ag,bg=engine(P5,G4,5),engine(Pb5,Q4,5)
            assert ag[1]==bg[1]==0
            G5,Q5=G4+(-ag[0],),Q4+(-bg[0],)
            S5=S4+[add(ns5,integer(-5*d))]
            ns6=base(P5,S5,6)
            e_box=box(ns6,6,41)
            if not e_box:
                continue
            for eta6,theta6,a6,b6 in update(A,B,engine(P5+[Z],G5,6),engine(Pb5+[Z],Q5,6),6):
                admissible=[(f,ea,eb) for f,ea,eb in options
                            if (a6==Z or r+s*(eb-f)==0) and (b6==Z or r+s*(ea-f)==0)]
                if not admissible:continue
                counters["row6"]+=1
                A6,B6=A+(a6,),B+(b6,)
                G6z,Q6z=G5+(eta6,),Q5+(theta6,)
                ref7=update(A6,B6,engine(P5+[Z,Z],G6z,7),engine(Pb5+[Z,Z],Q6z,7),7)
                for eta7,theta7,a7,b7 in ref7:
                    counters["row7"]+=1
                    A7,B7=A6+(a7,),B6+(b7,)
                    for f,ea,eb in admissible:
                        # Compute reverse third cofactor coefficients directly for rho7 compatibility.
                        eg=(ea-f)*r-s*f*(ea-f)
                        eq=(eb-f)*r-s*f*(eb-f)
                        choices=e_box
                        if b7!=Z:
                            choices=(eg//s,) if eg%s==0 and eg//s in e_box else ()
                        if a7!=Z:
                            choices=(eq//s,) if eq%s==0 and eq//s in choices else ()
                        for e in choices:
                            P8= P5+[integer(e),integer(f),I]
                            Pb8=Pb5+[integer(e),integer(f),I]
                            S6=S5+[add(ns6,integer(-6*e))]
                            s7=add(base(P8,S6,7),integer(-7*f))
                            if not house(s7,16385,32):
                                continue
                            s8=add(base(P8,S6+[s7],8),integer(-8))
                            if not house(s8,65537,64):
                                continue
                            counters["newton8_top_option"]+=1
                            G6,Q6=G5+(eta6-e,),Q5+(theta6-e,)
                            G7=G6+(eta7-f+e*(L-g1),)
                            Q7=Q6+(theta7-f+e*(L-q1),)
                            assert add(engine(P8,G6,7),integer(G7[7]))==a7
                            assert add(engine(Pb8,Q6,7),integer(Q7[7]))==b7
                            for _,_,a8,b8 in update(A7,B7,engine(P8,G7,8),engine(Pb8,Q7,8),8):
                                good=False
                                for xa,xb in TOP:
                                    if (B[3]!=Z and xa) or (A[3]!=Z and xb):continue
                                    reversed_P=P8[::-1]
                                    reversed_Pb=Pb8[::-1]
                                    def projection(poly,x1,x3):
                                        prefix=[1]
                                        for k,target in ((1,x1),(2,0),(3,x3)):
                                            y=engine(poly,prefix,k)
                                            assert y[1]==0
                                            prefix.append(target-y[0])
                                        return add(mulc(prefix[2],poly[6]),mulc(prefix[3],poly[5]))[1]
                                    if a8!=Z and projection(reversed_Pb,eb,xb):
                                        continue
                                    if b8!=Z and projection(reversed_P,ea,xa):
                                        continue
                                    good=True
                                    break
                                if good:
                                    counters["row8_top_options"]+=1
                                    key=(D,L,Ta,r,Tb,s,c,d,e,f)
                                    if key not in found:
                                        found[key]=env["annulus"](P8)
    return counters,found


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("front",type=Path)
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    raw=a.front.read_bytes()
    rows=[tuple(json.loads(line)) for line in raw.splitlines()]
    totals=Counter()
    all_rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for counts,found in pool.map(solve,rows,chunksize=8):
            totals.update(counts)
            for key,code in found.items():
                assert key not in all_rows or all_rows[key]==code
                all_rows[key]=code
    claimed=a.reference.read_bytes()
    actual=b"".join(json.dumps((row,code),separators=(",",":")).encode()+b"\n"
                    for row,code in sorted(all_rows.items()))
    assert actual==claimed
    expected=json.loads(a.reference.with_name("tail.json").read_bytes())
    assert dict(totals)==expected["counts"]
    out=dict(status="verified",counts=dict(totals),rows=len(all_rows),
             front_sha256=sha256(raw).hexdigest(),rows_sha256=sha256(claimed).hexdigest())
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__=="__main__":
    main()
