#!/usr/bin/env python3
"""Independent integral-basis S25 replay using complete integral-basis coordinates and direct common shifts."""

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

    def d_box(point_so_far):
        m,n=point_so_far
        if n:
            radius=(514-isqrt(D*(4*radicand_factor*n)**2)-1)//4
        else:
            radius=(514-1)//4
        u=2*m+t*n
        if radius<0:
            return range(0)
        return range(max(-120,(u-radius+7)//8),min(120,(u+radius)//8)+1)

    def root(P):
        ps=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))
        for x,y in ps:
            v=Z
            for j,c in enumerate(P):
                v=add(v,mulc(x**j*y**(8-j),c))
            for side,z in enumerate((v,conjugate(v))):
                if sign(z)<=0:
                    return "root_positive_reject",{"point":[x,y],"side":side,"value":z}
        H=[Z]*17
        for i,a in enumerate(P):
            for j,b in enumerate(map(conjugate,P)):
                H[i+j]=add(H[i+j],mul(a,b))
        assert all(x[1]==0 for x in H)
        norm=[x[0] for x in H]
        for side in ("outer","inner"):
            seq=norm if side=="outer" else norm[::-1]
            asc=[x<<j for j,x in enumerate(seq)]
            for n in range(16,0,-1):
                if abs(asc[0])>=abs(asc[-1]):
                    return f"root_{side}_reject",{"failure_degree":n,"norm":norm}
                low,high=asc[0],asc[-1]
                b=[high*asc[i+1]-low*asc[n-i-1] for i in range(n)]
                content=gcd(*b)
                asc=[x//content for x in b]
        return "root_pass",{"norm":norm}

    return locals()


TOP=((0,0),(1,0),(0,1))


def solve(D):
    e=basis(D)
    t, sign, house, choose, conj=e["t"],e["sign"],e["both_abs"],e["choices"],e["conjugate"]
    mul, update, engine, newton, cross=e["mul"],e["update"],e["engine"],e["newton_base"],e["cross"]
    counts=Counter()
    rows={}
    for r in range(1,9):
        for m in range(-13,14):
            a=(m,r)
            if not house(a,17,2):continue
            counts["a_cores"]+=1
            ab=conj(a)
            (g2,a2),=choose(a)
            (q2,b2),=choose(ab)
            assert add(a2,b2)==I
            A3,B3=(I,Z,a2,Z),(I,Z,b2,Z)
            for rho,sigma,a4,b4 in update(A3,B3,mulc(g2,a),mulc(q2,ab),4,(I,)):
                counts["a_row4"]+=1
                A4,B4=A3+(a4,),B3+(b4,)
                for c in range(-10,11):
                    G3,Q3=(1,0,g2,-c),(1,0,q2,-c)
                    S3=[integer(8),Z,mulc(-2,a),integer(-3*c)]
                    s4base=mulc(2,mul(a,a))
                    d_values=e["d_box"](s4base)
                    for s in range(-r,r+1):
                        if not s or r%s:continue
                        exceptional=s==r*c
                        for u in range(c*m-67,c*m+68):
                            b=(u,s)
                            bb=conj(b)
                            s5=mulc(5,add(mulc(c,a),mulc(-1,b)))
                            if not house(b,185,2) or not house(s5,1025,8):
                                continue
                            counts["direct_b_newton5"]+=1
                            (g3z,a3r),=choose(b)
                            (q3z,b3r),=choose(bb)
                            assert add(a3r,b3r)==I
                            for x,y in TOP:
                                for f in (0,1) if exceptional else (0,):
                                    if (a4!=Z and y!=f) or (b4!=Z and x!=f):continue
                                    for ee in range(-41,42) if exceptional else (0,):
                                        R3=[I,integer(f),integer(ee),b]
                                        RB3=[I,integer(f),integer(ee),bb]
                                        def reversed_g(poly,k1,kr3,ar3):
                                            U=[1]
                                            for k,z in ((1,integer(k1)),(2,Z),(3,ar3)):
                                                candidate=add(z,mulc(-1,engine(poly,U,k)))
                                                assert candidate[1]==0
                                                U.append(candidate[0])
                                            return tuple(U)
                                        Rg=reversed_g(R3,x,g3z,a3r)
                                        Rq=reversed_g(RB3,y,q3z,b3r)
                                        az=(r+s*Rq[3])!=0
                                        bz=(r+s*Rg[3])!=0
                                        if az and bz:continue
                                        counts["top_options"]+=1
                                        A3r,B3r=(I,integer(x),Z,a3r),(I,integer(y),Z,b3r)
                                        for d in d_values:
                                            if (az and d!=rho) or (bz and d!=sigma): continue
                                            counts["newton4_top_d"]+=1
                                            G4,Q4=G3+(rho-d,),Q3+(sigma-d,)
                                            P=[I,Z,a,integer(c),integer(d),b,integer(ee),integer(f),I]
                                            Pb=list(map(conj,P))
                                            S=S3+[add(s4base,integer(-4*d)),s5]
                                            for k,pair in ((6,(4097,16)),(7,(16385,32)),(8,(65537,64))):
                                                sk=add(newton(P,S,k),mulc(-k,P[k]))
                                                if not house(sk,*pair):break
                                                S.append(sk)
                                            if len(S)!=9:continue
                                            counts["newton8_top_d"]+=1
                                            for g5,q5,a5,b5 in update(A4,B4,engine(P,G4,5),engine(Pb,Q4,5),5):
                                                G5,Q5=G4+(g5,),Q4+(q5,)
                                                A5,B5=A4+(a5,),B4+(b5,)
                                                for g6,q6,a6,b6 in update(A5,B5,engine(P,G5,6),engine(Pb,Q5,6),6,(I,)):
                                                    if (az and a6!=Z) or (bz and b6!=Z):continue
                                                    counts["row6"]+=1
                                                    G6,Q6=G5+(g6,),Q5+(q6,)
                                                    A6,B6=A5+(a6,),B5+(b6,)
                                                    for g7,q7,a7,b7 in update(A6,B6,engine(P,G6,7),engine(Pb,Q6,7),7):
                                                        counts["row7"]+=1
                                                        G7,Q7=G6+(g7,),Q6+(q7,)
                                                        A7,B7=A6+(a7,),B6+(b7,)
                                                        for _,_,a8,b8 in update(A7,B7,engine(P,G7,8),engine(Pb,Q7,8),8):
                                                            counts["row8"]+=1
                                                            A8,B8=A7+(a8,),B7+(b8,)
                                                            active=[(Rg,Rq,A3r,B3r)]
                                                            for k in range(4,9):
                                                                nxt=[]
                                                                for U,V,C,E in active:
                                                                    R,RS=P[::-1],Pb[::-1]
                                                                    for ui,vi,ak,bk in update(C,E,engine(R,U,k),engine(RS,V,k),k):
                                                                        if (A8[k]!=Z and bk!=Z) or (B8[k]!=Z and ak!=Z):
                                                                            continue
                                                                        nxt.append((U+(ui,),V+(vi,),C+(ak,),E+(bk,)))
                                                                active=nxt
                                                                if not active:break
                                                            if not active:continue
                                                            counts["coupled_row8"]+=1
                                                            mother=(D,2*m+t*r,r,c,d,2*u+t*s,s,ee,f)
                                                            if mother not in rows:
                                                                rows[mother]=e["root"](P)
                                                                counts["mothers_exception" if exceptional else "mothers_generic"]+=1
    return counts, rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,289) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    totals=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for counts,local in pool.map(solve,domains,chunksize=1):
            totals.update(counts)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
    reference=a.reference.read_bytes()
    assert {tuple(row):code for row,code in map(json.loads,reference.splitlines())}=={
        row:code for row,(code,witness) in rows.items()}
    other=json.loads(a.reference.with_name("entry.json").read_bytes())["counts"]
    assert all(v==other[k] for k,v in totals.items() if k!="direct_b_newton5")
    out=dict(status="verified",totals=dict(totals),rows=len(rows),
             reference_sha256=sha256(reference).hexdigest(),
             witness=[(row,*data) for row,data in sorted(rows.items())])
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="witness"},sort_keys=True))


if __name__=="__main__":
    main()
