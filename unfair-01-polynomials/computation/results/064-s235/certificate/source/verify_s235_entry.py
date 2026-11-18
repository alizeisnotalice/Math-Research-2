#!/usr/bin/env python3
"""Independent integral-basis S235 replay using complete integral-basis coordinates and direct high-coordinate tests."""

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


def high_options(env,a,r,L,u,a4,b4):
    D,eps,tau=env["D"],env["t"],env["radicand_factor"]
    shifts,house,conj,engine=env["choices"],env["both_abs"],env["conjugate"],env["engine"]
    for irr in range(-r,r+1):
        if not irr or r%irr or irr==r*u:continue
        if (40*tau*(r*u-irr))**2*D>=2050**2:continue
        e=L*r//irr
        if not -41<=e<=41:continue
        base=(0,irr)
        b=conj(base)
        candidates={r//irr+shifts(base)[0][0],r//irr+shifts(b)[0][0]}
        for v in sorted(candidates):
            d=add(base,integer(v))
            if not house(d,185,2):continue
            for x,y in TOP:
                if (a4!=Z and y) or (b4!=Z and x):continue
                yield d,irr,e,0,x,y,False
    if not u: return
    irr=r*u
    if (tau*irr)**2*D>=185**2:return
    d0=mulc(u,a)
    d0b=conj(d0)
    rA,rB=shifts(d0)[0][1],shifts(d0b)[0][1]
    for x,y in TOP:
        for f in (0,1):
            if (a4!=Z and y!=f) or (b4!=Z and x!=f):continue
            for e in range(-41,42):
                P0,Q0=[I,integer(f),integer(e),d0],[I,integer(f),integer(e),d0b]
                U=[1]
                V=[1]
                for k,(ak,bk) in enumerate(((integer(x),integer(y)),(Z,Z),(rA,rB)),1):
                    one=add(ak,mulc(-1,engine(P0,U,k)))
                    two=add(bk,mulc(-1,engine(Q0,V,k)))
                    assert one[1]==two[1]==0
                    U.append(one[0])
                    V.append(two[0])
                cases=set()
                for rev in (U,V):
                    q=r+L*r*rev[1]+irr*rev[3]
                    if q%irr==0:
                        cases.add(q//irr)
                for v in sorted(cases):
                    d=add(d0,integer(v))
                    if house(d,185,2):yield d,irr,e,f,x,y,True


def solve(D):
    E=basis(D)
    eps=E["t"]
    choices,conj,house=E["choices"],E["conjugate"],E["both_abs"]
    engine,update,newton=E["engine"],E["update"],E["newton_base"]
    rows={}
    counts=Counter()
    for L in (-4,-3,-2,-1,1,2,3,4):
        for r in range(1,9):
            for m in range(-21,22):
                a=(m,r)
                S2=add(integer(L*L),mulc(-2,a))
                if not house(S2,17,1) or not house(a,83,2):continue
                counts["a_cores"]+=1
                ab=conj(a)
                (rho,a2),=choices(a)
                (sigma,b2),=choices(ab)
                assert add(a2,b2)==I
                A3,B3=(I,Z,a2,Z),(I,Z,b2,Z)
                G2,Q2=(1,-L,L*L+rho),(1,-L,L*L+sigma)
                for u in range((-L**3)//3-12,(-L**3)//3+13):
                    s3=integer(-L**3-3*u)
                    if not house(s3,65,2):continue
                    b=add(mulc(L,a),integer(u))
                    if not house(b,185,2):continue
                    counts["newton3_b"]+=1
                    bb=conj(b)
                    P3,Pb3=[I,integer(L),a,b],[I,integer(L),ab,bb]
                    G3,Q3=G2+(-u-L*G2[2],),Q2+(-u-L*Q2[2],)
                    S3=[integer(8),integer(-L),S2,s3]
                    ns4=newton(P3,S3,4)
                    c_values=E["d_box"](ns4)
                    for rho4,sigma4,a4,b4 in update(A3,B3,engine(P3+[Z],G3,4),
                                                    engine(Pb3+[Z],Q3,4),4,(I,)):
                        counts["row4"]+=1
                        A4,B4=A3+(a4,),B3+(b4,)
                        for d,irr,e,f,x,y,exceptional in high_options(E,a,r,L,u,a4,b4):
                            counts["high_options"]+=1
                            dbl=conj(d)
                            R3,Rb3=[I,integer(f),integer(e),d],[I,integer(f),integer(e),dbl]
                            (ks,az),=choices(d)
                            (ls,bz),=choices(dbl)
                            X3,Y3=(I,integer(x),Z,az),(I,integer(y),Z,bz)
                            U,V=[1],[1]
                            for k,(ak,bk) in enumerate(zip(X3[1:],Y3[1:]),1):
                                one=add(ak,mulc(-1,engine(R3,U,k)))
                                two=add(bk,mulc(-1,engine(Rb3,V,k)))
                                assert one[1]==two[1]==0
                                U.append(one[0])
                                V.append(two[0])
                            a_zero=r+L*r*V[1]+irr*V[3]!=0
                            b_zero=r+L*r*U[1]+irr*U[3]!=0
                            assert not(a_zero and b_zero)
                            for c in c_values:
                                G4,Q4=G3+(rho4-c,),Q3+(sigma4-c,)
                                if a_zero and r*G4[4]+L*r*G4[3]-L*irr:continue
                                if b_zero and r*Q4[4]+L*r*Q4[3]-L*irr:continue
                                P=P3+[integer(c),d,integer(e),integer(f),I]
                                Pb=Pb3+[integer(c),dbl,integer(e),integer(f),I]
                                S=S3+[add(ns4,integer(-4*c))]
                                for k,bound in ((5,(1025,8)),(6,(4097,16)),(7,(16385,32)),(8,(65537,64))):
                                    power=add(newton(P,S,k),mulc(-k,P[k]))
                                    if not house(power,*bound):break
                                    S.append(power)
                                if len(S)!=9:continue
                                counts["newton8"]+=1
                                for g5,q5,a5,b5 in update(A4,B4,engine(P,G4,5),engine(Pb,Q4,5),5):
                                    G5,Q5=G4+(g5,),Q4+(q5,)
                                    A5,B5=A4+(a5,),B4+(b5,)
                                    for g6,q6,a6,b6 in update(A5,B5,engine(P,G5,6),engine(Pb,Q5,6),6,(I,)):
                                        if (a_zero and a6!=Z) or (b_zero and b6!=Z):continue
                                        counts["row6"]+=1
                                        G6,Q6=G5+(g6,),Q5+(q6,)
                                        A6,B6=A5+(a6,),B5+(b6,)
                                        for g7,q7,a7,b7 in update(A6,B6,engine(P,G6,7),engine(Pb,Q6,7),7):
                                            G7,Q7=G6+(g7,),Q6+(q7,)
                                            A7,B7=A6+(a7,),B6+(b7,)
                                            for _,_,a8,b8 in update(A7,B7,engine(P,G7,8),engine(Pb,Q7,8),8):
                                                counts["row8"]+=1
                                                A8,B8=A7+(a8,),B7+(b8,)
                                                R,Rb=P[::-1],Pb[::-1]
                                                live=((tuple(U),tuple(V),X3,Y3),)
                                                for k in range(4,9):
                                                    nxt=[]
                                                    for X,Y,W,Z_ in live:
                                                        for ui,vi,ak,bk in update(W,Z_,engine(R,X,k),engine(Rb,Y,k),k):
                                                            if (A8[k]!=Z and bk!=Z) or (B8[k]!=Z and ak!=Z):
                                                                continue
                                                            nxt.append((X+(ui,),Y+(vi,),W+(ak,),Z_+(bk,)))
                                                    live=nxt
                                                    if not live:break
                                                if not live:continue
                                                counts["coupled_row8"]+=1
                                                row=(D,L,2*m+eps*r,r,u,c,2*d[0]+eps*irr,irr,e,f)
                                                if row not in rows:
                                                    rows[row]=E["root"](P)
                                                    counts["mothers_exception" if exceptional else "mothers_generic"]+=1
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
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for c,rs in pool.map(solve,domains,chunksize=1):
            totals.update(c)
            assert rows.keys().isdisjoint(rs)
            rows.update(rs)
    raw=a.reference.read_bytes()
    actual=b"".join(json.dumps((row,code),separators=(",",":")).encode()+b"\n"
                    for row,(code,_) in sorted(rows.items()))
    assert actual==raw
    expected=json.loads(a.reference.with_name("entry.json").read_bytes())["counts"]
    assert all(expected[k]==v for k,v in totals.items())
    out=dict(status="verified",counts=dict(totals),mothers=len(rows),
             rows_sha256=sha256(actual).hexdigest(),
             witness=[(row,code,detail) for row,(code,detail) in sorted(rows.items())])
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="witness"},sort_keys=True))


if __name__=="__main__":
    main()
