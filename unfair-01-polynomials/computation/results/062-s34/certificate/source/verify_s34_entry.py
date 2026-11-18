#!/usr/bin/env python3
"""Independent integral-basis S34 replay using complete integral-basis coordinates and direct common shifts."""

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

    def c_range(U):
        m,n=U
        R=(2050-isqrt(D*(8*radicand_factor*n)**2)-1)//8
        if R<0:return range(0)
        x=2*m+t*n
        return range(max(-92,(x-R+9)//10),min(92,(x+R)//10)+1)

    def check_coupling(P,G,Q,A,B,forward_a,forward_b,max_row):
        pb=list(map(conjugate,P))
        live=((G,Q,A,B),)
        for k in range(4,max_row+1):
            candidates=[]
            for U,V,C,E in live:
                for ui,vi,ck,ek in update(C,E,engine(P,U,k),engine(pb,V,k),k):
                    if (forward_b[k]!=Z and ck!=Z) or (forward_a[k]!=Z and ek!=Z):
                        continue
                    candidates.append((U+(ui,),V+(vi,),C+(ck,),E+(ek,)))
            if not candidates:
                return False
            live=candidates
        return True

    return locals()


TOP=((0,0),(1,0),(0,1))


def solve(D):
    en=basis(D)
    t=en["t"]
    mul,conj,shifts,house=en["mul"],en["conjugate"],en["choices"],en["both_abs"]
    update,engine,cross,newton=en["update"],en["engine"],en["cross"],en["newton_base"]
    coupled,crange,root=en["check_coupling"],en["c_range"],en["root"]
    counts=Counter()
    cases={}
    for L in (-4,-3,-2,-1,1,2,3,4):
        for M in range(-7,9) if abs(L)==1 else (0,):
            center=(-L**3+3*L*M)//3
            for r in range(1,11):
                for m in range(center-20,center+21):
                    a=(m,r)
                    a3s=add(integer(-L**3+3*L*M),mulc(-3,a))
                    if not house(a,185,2) or not house(a3s,65,2): continue
                    counts["a_cores"]+=1
                    ab=conj(a)
                    (rho,a3),=shifts(a)
                    (sigma,b3),=shifts(ab)
                    assert add(a3,b3)==I
                    S3=[integer(8),integer(-L),integer(L*L-2*M),a3s]
                    P3=[I,integer(L),integer(M),a]
                    b0=mulc(L,a)
                    P0=P3+[b0,Z,Z,Z,I]
                    Pb0=list(map(conj,P0))
                    ns4=newton(P0,S3,4)
                    s4z=add(ns4,mulc(-4,b0))
                    assert s4z[1]==0
                    middle=s4z[0]//4
                    s5z=newton(P0,S3+[s4z],5)
                    if (8*en["radicand_factor"]*s5z[1])**2*D>=2050**2:
                        continue
                    s6z=newton(P0,S3+[s4z,s5z],6)
                    if (16*en["radicand_factor"]*s6z[1])**2*D>=8194**2:
                        continue
                    for a2,b2 in TOP if abs(L)==1 else ((0,0),):
                        g2=L*L-M+a2
                        q2=L*L-M+b2
                        g3=rho+M*L-L*g2
                        q3=sigma+M*L-L*q2
                        G3,Q3=(1,-L,g2,g3),(1,-L,q2,q3)
                        g4,q4=engine(P0,G3,4),engine(Pb0,Q3,4)
                        assert g4[1]==q4[1]==0
                        G4,Q4=G3+(-g4[0],),Q3+(-q4[0],)
                        A4,B4=(I,Z,integer(a2),a3,Z),(I,Z,integer(b2),b3,Z)
                        for g5,q5,a5,b5 in update(A4,B4,engine(P0,G4,5),engine(Pb0,Q4,5),5):
                            A5,B5=A4+(a5,),B4+(b5,)
                            G5,Q5=G4+(g5,),Q4+(q5,)
                            for _,_,a6,b6 in update(A5,B5,engine(P0,G5,6),engine(Pb0,Q5,6),6,(I,)):
                                A6,B6=A5+(a6,),B5+(b6,)
                                for x,y in TOP:
                                    for p,q in TOP:
                                        if a2*q or b2*p:continue
                                        for ee in range(-9,10):
                                            U1,V1=x-ee,y-ee
                                            if (b5!=Z and 1+L*U1) or (a5!=Z and 1+L*V1):
                                                continue
                                            U2z,V2z=p-ee*U1,q-ee*V1
                                            w=(V1+L*V2z) if a6!=Z else (U1+L*U2z)
                                            if w%L:continue
                                            d=w//L
                                            if not -41<=d<=41:continue
                                            U2,V2=U2z-d,V2z-d
                                            if (b6!=Z and U1+L*U2) or (a6!=Z and V1+L*V2):
                                                continue
                                            Rbase=[I,integer(ee),integer(d),Z,b0,a,integer(M),integer(L),I]
                                            Gtop0,Qtop0=(1,U1,U2,-d*U1-ee*U2),(1,V1,V2,-d*V1-ee*V2)
                                            Atop,Btop=(I,integer(x),integer(p),Z),(I,integer(y),integer(q),Z)
                                            if not coupled(Rbase,Gtop0,Qtop0,Atop,Btop,A6,B6,6):
                                                continue
                                            counts["early_reverse_top"]+=1
                                            for u in range(middle-17,middle+18):
                                                b=add(b0,integer(u))
                                                if not house(b,1921,16):continue
                                                s4=add(s4z,integer(-4*u))
                                                if not house(s4,257,4):continue
                                                P4=P3+[b]
                                                S4=S3+[s4]
                                                ns5=newton(P4,S4,5)
                                                dummy=P4+[Z,integer(d),integer(ee),I]
                                                ns6=add(newton(dummy,S4+[ns5],6),integer(-6*d))
                                                ns7=add(newton(dummy,S4+[ns5,ns6],7),integer(-7*ee))
                                                if (32*en["radicand_factor"]*ns7[1])**2*D>=32770**2:
                                                    continue
                                                c_range=crange(ns5)
                                                if u+M*g2 or u+M*q2:
                                                    W=(V2,V1) if u+M*g2 else (U2,U1)
                                                    if W[0]%L:continue
                                                    candidate=W[0]//L-d*W[1]-ee*W[0]
                                                    c_range=(candidate,) if candidate in c_range else ()
                                                for c in c_range:
                                                    U3,V3=-c-d*U1-ee*U2,-c-d*V1-ee*V2
                                                    if (u+M*g2 and V2+L*V3) or (u+M*q2 and U2+L*U3):
                                                        continue
                                                    P=P4+[integer(c),integer(d),integer(ee),I]
                                                    Pb=list(map(conj,P))
                                                    S=S4+[add(ns5,integer(-5*c))]
                                                    for k,bound in ((6,(4097,16)),(7,(16385,32)),(8,(65537,64))):
                                                        v=add(newton(P,S,k),mulc(-k,P[k]))
                                                        if not house(v,*bound):break
                                                        S.append(v)
                                                    if len(S)!=9:continue
                                                    G4r,Q4r=G3+(G4[4]-u,),Q3+(Q4[4]-u,)
                                                    H5,J5=add(a5,mulc(-1,engine(P,G4r,5))),add(b5,mulc(-1,engine(Pb,Q4r,5)))
                                                    assert H5[1]==J5[1]==0
                                                    G5r,Q5r=G4r+(H5[0],),Q4r+(J5[0],)
                                                    H6,J6=add(a6,mulc(-1,engine(P,G5r,6))),add(b6,mulc(-1,engine(Pb,Q5r,6)))
                                                    assert H6[1]==J6[1]==0
                                                    G6,Q6=G5r+(H6[0],),Q5r+(J6[0],)
                                                    for g7,q7,a7,b7 in update(A6,B6,engine(P,G6,7),engine(Pb,Q6,7),7):
                                                        G7,Q7=G6+(g7,),Q6+(q7,)
                                                        A7,B7=A6+(a7,),B6+(b7,)
                                                        for _,_,a8,b8 in update(A7,B7,engine(P,G7,8),engine(Pb,Q7,8),8):
                                                            if not coupled(P[::-1],(1,U1,U2,U3),(1,V1,V2,V3),
                                                                           Atop,Btop,A7+(a8,),B7+(b8,),8):
                                                                continue
                                                            counts["coupled_row8"]+=1
                                                            row=(D,L,M,2*m+t*r,r,u,c,d,ee)
                                                            if row not in cases:
                                                                cases[row]=root(P)[0]
    return counts,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,470) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    totals=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for ct,local in pool.map(solve,domains,chunksize=1):
            totals.update(ct)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
    payload=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in sorted(rows.items()))
    reference=a.reference.read_bytes()
    assert payload==reference
    baseline=json.loads(a.reference.with_name("entry.json").read_bytes())["counts"]
    assert all(totals[k]==baseline[k] for k in ("a_cores","coupled_row8"))
    out=dict(status="verified",counts=dict(totals),mothers=len(rows),rows_sha256=sha256(payload).hexdigest(),
             codes=dict(Counter(rows.values())))
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__=="__main__":
    main()
