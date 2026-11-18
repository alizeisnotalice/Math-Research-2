#!/usr/bin/env python3
"""Independent integral-basis S245 replay using complete integral-basis coordinates and direct common shifts."""

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


def high_options(E,a,r,b,t,j,a4,b4):
    D,eps,tau=E["D"],E["t"],E["radicand_factor"]
    shifts,house,conj,engine=E["choices"],E["both_abs"],E["conjugate"],E["engine"]
    exceptional=t==r*b
    if not exceptional:
        if not -9<=j<=9:return
        if (40*tau*(r*b-t))**2*D>=2050**2:return
        db=(0,t)
        candidates={r//t+shifts(db)[0][0],r//t+shifts(conj(db))[0][0]}
        for v in sorted(candidates):
            d=add(db,integer(v))
            if not house(d,185,2):continue
            for x,y in TOP:
                if (a4!=Z and y) or (b4!=Z and x):continue
                yield d,0,j,x,y,False
        return
    assert abs(b)==1
    d0=mulc(b,a)
    d0b=conj(d0)
    one,two=shifts(d0)[0][1],shifts(d0b)[0][1]
    for f in (j,j+1):
        if not -9<=f<=9:continue
        for x,y in TOP:
            if (a4!=Z and y-f+j) or (b4!=Z and x-f+j):continue
            for e in range(-41,42):
                R,Rb=[I,integer(f),integer(e),d0],[I,integer(f),integer(e),d0b]
                U,V=[1],[1]
                for k,(ak,bk) in enumerate(((integer(x),integer(y)),(Z,Z),(one,two)),1):
                    newA=add(ak,mulc(-1,engine(R,U,k)))
                    newB=add(bk,mulc(-1,engine(Rb,V,k)))
                    assert newA[1]==newB[1]==0
                    U.append(newA[0])
                    V.append(newB[0])
                trial={r//t+j*U[2]+U[3],r//t+j*V[2]+V[3]}
                for v in sorted(trial):
                    if -25<=v<=25:
                        d=add(d0,integer(v))
                        if house(d,185,2):yield d,e,f,x,y,True


def solve(D):
    env=basis(D)
    eps=env["t"]
    choose,conj,house=env["choices"],env["conjugate"],env["both_abs"]
    engine,update,newton=env["engine"],env["update"],env["newton_base"]
    mul=env["mul"]
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
            a_sq=mul(a,a)
            for t in range(-r,r+1):
                if not t or r%t:continue
                for j in range(-10,10):
                    if not j:continue
                    s=t*j
                    c0=(0,s)
                    cb0=conj(c0)
                    base=add(mulc(2,a_sq),mulc(-4,c0))
                    tr=2*base[0]+eps*base[1]
                    lo=-((128-tr)//8)
                    hi=(tr+128)//8
                    cs=[]
                    for v in range(lo,hi+1):
                        c=(v,s)
                        s4=add(base,integer(-4*v))
                        if house(c,1921,16) and house(s4,257,4):
                            cs.append((v,c,s4))
                    if not cs:continue
                    counts["core_a_t_j"]+=1
                    for rho,sigma,a4,b4 in update(A3,B3,add(c0,mulc(g2,a)),
                                                  add(cb0,mulc(q2,ab)),4,(I,)):
                        counts["row4_normalized"]+=1
                        A4,B4=A3+(a4,),B3+(b4,)
                        for b in range(-10,11):
                            P3,Pb3=[I,Z,a,integer(b)],[I,Z,ab,integer(b)]
                            G3,Q3=(1,0,g2,-b),(1,0,q2,-b)
                            S3=[integer(8),Z,mulc(-2,a),integer(-3*b)]
                            for d,e,f,x,y,exceptional in high_options(env,a,r,b,t,j,a4,b4):
                                db=conj(d)
                                s5=mulc(5,add(mulc(b,a),mulc(-1,d)))
                                if not house(s5,1025,8):
                                    continue
                                counts["high_options"]+=1
                                R3,Rb3=[I,integer(f),integer(e),d],[I,integer(f),integer(e),db]
                                (g3r,ar3),=choose(d)
                                (q3r,br3),=choose(db)
                                X3,Y3=(I,integer(x),Z,ar3),(I,integer(y),Z,br3)
                                U,V=[1],[1]
                                for k,(ak,bk) in enumerate(zip(X3[1:],Y3[1:]),1):
                                    nextA=add(ak,mulc(-1,engine(R3,U,k)))
                                    nextB=add(bk,mulc(-1,engine(Rb3,V,k)))
                                    assert nextA[1]==nextB[1]==0
                                    U.append(nextA[0])
                                    V.append(nextB[0])
                                revA=r+s*U[2]+t*U[3]
                                revB=r+s*V[2]+t*V[3]
                                assert not(revA and revB)
                                for v,c,s4 in cs:
                                    G4,Q4=G3+(rho-v,),Q3+(sigma-v,)
                                    if (revB and r*G4[4]+s*g2) or (revA and r*Q4[4]+s*q2):continue
                                    P=P3+[c,d,integer(e),integer(f),I]
                                    Pb=Pb3+[conj(c),db,integer(e),integer(f),I]
                                    S=S3+[s4,s5]
                                    for k,bound in ((6,(4097,16)),(7,(16385,32)),(8,(65537,64))):
                                        q=add(newton(P,S,k),mulc(-k,P[k]))
                                        if not house(q,*bound):break
                                        S.append(q)
                                    if len(S)!=9:continue
                                    counts["newton8"]+=1
                                    for g5,q5,a5,b5 in update(A4,B4,engine(P,G4,5),engine(Pb,Q4,5),5):
                                        A5,B5=A4+(a5,),B4+(b5,)
                                        G5,Q5=G4+(g5,),Q4+(q5,)
                                        for g6,q6,a6,b6 in update(A5,B5,engine(P,G5,6),engine(Pb,Q5,6),6,(I,)):
                                            if (revB and a6!=Z) or (revA and b6!=Z):continue
                                            counts["row6"]+=1
                                            A6,B6=A5+(a6,),B5+(b6,)
                                            G6,Q6=G5+(g6,),Q5+(q6,)
                                            for g7,q7,a7,b7 in update(A6,B6,engine(P,G6,7),engine(Pb,Q6,7),7):
                                                A7,B7=A6+(a7,),B6+(b7,)
                                                G7,Q7=G6+(g7,),Q6+(q7,)
                                                for _,_,a8,b8 in update(A7,B7,engine(P,G7,8),engine(Pb,Q7,8),8):
                                                    counts["row8"]+=1
                                                    A8,B8=A7+(a8,),B7+(b8,)
                                                    R,Rb=P[::-1],Pb[::-1]
                                                    live=((tuple(U),tuple(V),X3,Y3),)
                                                    for k in range(4,9):
                                                        nxt=[]
                                                        for X,Y,W,Z_ in live:
                                                            for a_,b_,ak,bk in update(W,Z_,engine(R,X,k),engine(Rb,Y,k),k):
                                                                if (A8[k]!=Z and bk!=Z) or (B8[k]!=Z and ak!=Z):
                                                                    continue
                                                                nxt.append((X+(a_,),Y+(b_,),W+(ak,),Z_+(bk,)))
                                                        live=nxt
                                                        if not live:break
                                                    if not live:continue
                                                    counts["coupled_row8"]+=1
                                                    row=(D,2*m+eps*r,r,b,2*v+eps*s,s,2*d[0]+eps*t,t,e,f)
                                                    if row not in rows:
                                                        rows[row]=env["root"](P)
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
