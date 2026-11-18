#!/usr/bin/env python3
"""Independent integral-basis S2345 replay using complete integral-basis coordinates and direct common shifts."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import gcd, isqrt
from functools import lru_cache
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
        for i in range(max(1,k-len(G)+1),min(k,len(P)-1)+1):
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

    return locals()

TOP=((0,0),(1,0),(0,1))


@lru_cache
def geometry(D,L,r):
    """Enumerate integer source boxes directly from f, and then derive j."""
    tau=1 if D%4==1 else 2
    all_B=range((-32-L**3+2)//3,(32-L**3)//3+1)
    bound=isqrt((65*65-1)//(9*tau*tau*D))
    result={}
    for f in range(-9,10):
        es=tuple(e for e in range(-41,42) if abs(f*f-2*e)<=16)
        for j in (f,f-1):
            if not j:continue
            for t in range(-bound,bound+1):
                if not t:continue
                exceptional_B=(t*(1-L*j))//r
                if r*exceptional_B==t*(1-L*j) and exceptional_B in all_B:
                    key=(exceptional_B,t,j)
                    result.setdefault(key,set()).update((e,f,True) for e in es)
                if j!=f or r%t or L*r//t not in es:continue
                for B in all_B:
                    if r*B!=t*(1-L*j):
                        result.setdefault((B,t,j),set()).add((L*r//t,f,False))
    by_B={}
    for (B,t,j),options in result.items():
        by_B.setdefault(B,[]).append((t,j,sorted(options)))
    return by_B


def slope_range(E,base,slope,k,lo,hi):
    """Integer trace slab, after reserving both embedding differences."""
    numerator,denominator=BOUND[k]
    q=E["radicand_factor"]*denominator*base[1]
    radius=(2*numerator-isqrt(E["D"]*q*q)-1)//denominator
    if radius<0:return 1,0
    tr=2*base[0]+E["t"]*base[1]
    if slope<0:tr,slope=-tr,-slope
    if not slope:return (lo,hi) if abs(tr)<=radius else (1,0)
    step=2*slope
    return max(lo,-((radius+tr)//step)),min(hi,(radius-tr)//step)


def moments(P,E,last):
    S=[integer(8)]
    for k in range(1,last+1):
        S.append(add(E["newton_base"](P,S,k),mulc(-k,P[k])))
    return S


def high_initial(E,R,x):
    """Direct inverse convolution from the high-end coefficient targets."""
    shifts,engine=E["choices"],E["engine"]
    (g,third),=shifts(R[3])
    A=(I,integer(x),Z,third)
    U=[1]
    for k in range(1,4):
        new=add(A[k],mulc(-1,engine(R,U,k)))
        assert new[1]==0
        U.append(new[0])
    return tuple(U),A


def compatible(E,P,U,V,X,Y,A,B):
    rev=P[::-1]
    other=list(map(E["conjugate"],rev))
    states=[(U,V,X,Y)]
    for k in range(4,9):
        nxt=[]
        for G,Q,C,F in states:
            for g,h,ak,bk in E["update"](C,F,E["engine"](rev,G,k),E["engine"](other,Q,k),k):
                if (A[k]!=Z and bk!=Z) or (B[k]!=Z and ak!=Z):continue
                nxt.append((G+(g,),Q+(h,),C+(ak,),F+(bk,)))
        if not nxt:return False
        states=nxt
    return True


def solve(D):
    E=basis(D)
    eps=E["t"]
    choose,conj,house=E["choices"],E["conjugate"],E["both_abs"]
    engine,update,newton=E["engine"],E["update"],E["newton_base"]
    counts=Counter()
    rows={}
    for L in range(-9,10):
        if not L:continue
        for r in range(1,9):
            for m in range(-54,55):
                a=(m,r)
                if not house(add(mulc(2,a),integer(-L*L)),17,1):continue
                if not house(a,83,2):continue
                counts["a_cores"]+=1
                ab=conj(a)
                (pa,a2),=choose(a)
                (pb,b2),=choose(ab)
                assert add(a2,b2)==I
                A3,B3=(I,Z,a2,Z),(I,Z,b2,Z)
                for B,combos in geometry(D,L,r).items():
                    b=add(mulc(L,a),integer(B))
                    if not house(b,185,2):continue
                    P3=(I,integer(L),a,b)
                    H3=tuple(map(conj,P3))
                    g2,q2=L*L+pa,L*L+pb
                    G3,Q3=(1,-L,g2,-B-L*g2),(1,-L,q2,-B-L*q2)
                    S3=moments(P3,E,3)
                    if not house(S3[3],*BOUND[3]):continue
                    for t,j,tails in combos:
                        s=t*j
                        c0=(0,s)
                        s4=add(newton(P3,S3,4),mulc(-4,c0))
                        lo,hi=slope_range(E,s4,-4,4,-174,174)
                        if lo>hi:continue
                        possible=[]
                        for v in range(lo,hi+1):
                            if house((v,s),1921,16):
                                possible.append(v)
                        if not possible:continue
                        counts["core_LaBtj"]+=1
                        lo,hi=possible[0],possible[-1]
                        for rho,sigma,a4,b4 in update(A3,B3,engine(P3+(c0,),G3,4),
                                                      engine(H3+(conj(c0),),Q3,4),4,(I,)):
                            counts["row4_normalized"]+=1
                            A4,B4=A3+(a4,),B3+(b4,)
                            d0=(0,t)
                            for e,f,exceptional in tails:
                                R6=(I,integer(f),integer(e),d0,c0,b,a)
                                R6b=tuple(map(conj,R6))
                                for x,y in TOP:
                                    U0,X0=high_initial(E,R6,x)
                                    V0,Y0=high_initial(E,R6b,y)
                                    if ((a4!=Z and engine(R6b,V0,4)[1]) or
                                        (b4!=Z and engine(R6,U0,4)[1])):
                                        continue
                                    hA,hB=engine(R6,U0,6)[1],engine(R6b,V0,6)[1]
                                    translations={h//slope for h,slope in ((hA,t),(hB,-t)) if h%slope==0}
                                    for w in sorted(translations):
                                        d=(w,t)
                                        R=(I,integer(f),integer(e),d,c0,b,a,integer(L),I)
                                        SR=moments(R,E,3)
                                        if not house(SR[3],*BOUND[3]) or not house(d,185,2):continue
                                        for k in (4,5):
                                            SR.append(add(newton(R,SR,k),mulc(-k,R[k])))
                                        counts["high_options"]+=1
                                        Rb=tuple(map(conj,R))
                                        U,X=high_initial(E,R,x)
                                        V,Y=high_initial(E,Rb,y)
                                        assert X==X0 and Y==Y0 and add(X[3],Y[3])==I
                                        revA,revB=engine(R,U,6)[1],engine(Rb,V,6)[1]
                                        assert not (revA and revB)
                                        P0=P3+(c0,d,integer(e),integer(f),I)
                                        s5=moments(P0,E,5)[5]
                                        ll,hh=slope_range(E,s5,5*L,5,lo,hi)
                                        ll,hh=slope_range(E,SR[4],-4,4,ll,hh)
                                        ll,hh=slope_range(E,SR[5],5*f,5,ll,hh)
                                        for v in range(ll,hh+1):
                                            P=P3+((v,s),d,integer(e),integer(f),I)
                                            Pb=tuple(map(conj,P))
                                            G4,Q4=G3+(rho-v,),Q3+(sigma-v,)
                                            if revB and engine(P,G4+(0,),6)[1]:continue
                                            if revA and engine(Pb,Q4+(0,),6)[1]:continue
                                            if any(not house(z,*BOUND[k]) for poly in (P,P[::-1])
                                                   for k,z in enumerate(moments(poly,E,8)) if k>=2):
                                                continue
                                            counts["newton8_both"]+=1
                                            states=[(G4,Q4,A4,B4)]
                                            for k in range(5,9):
                                                nxt=[]
                                                for G,Q,A,F in states:
                                                    for g,h,ak,bk in update(A,F,engine(P,G,k),engine(Pb,Q,k),k):
                                                        if k==6 and ((revB and ak!=Z) or (revA and bk!=Z)):
                                                            continue
                                                        nxt.append((G+(g,),Q+(h,),A+(ak,),F+(bk,)))
                                                states=nxt
                                                if not states:break
                                                counts[f"row{k}"]+=len(states)
                                            for _,_,A,F in states:
                                                if not compatible(E,P,U,V,X,Y,A,F):continue
                                                counts["coupled_row8"]+=1
                                                row=(D,L,2*m+eps*r,r,B,2*v+eps*s,s,2*w+eps*t,t,e,f)
                                                if row not in rows:
                                                    rows[row]=E["root"](P)
                                                    counts["mothers_exception" if exceptional else "mothers_generic"]+=1
    return counts,rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=8)
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
    if actual!=raw:
        raise ValueError(dict(actual=len(rows),expected=len(raw.splitlines()),counts=totals,
                              missing=[x.decode() for x in sorted(set(raw.splitlines())-set(actual.splitlines()))[:10]],
                              extra=[x.decode() for x in sorted(set(actual.splitlines())-set(raw.splitlines()))[:10]]))
    expected=json.loads(a.reference.with_name("entry.json").read_bytes())["counts"]
    assert all(expected[k]==v for k,v in totals.items()), (totals,expected)
    details=[[row,code,data] for row,(code,data) in sorted(rows.items())]
    result=dict(status="verified",counts=totals,rows_sha256=sha256(actual).hexdigest(),details=details)
    a.out.write_text(json.dumps(result,sort_keys=True,separators=(",",":"))+"\n")
    print(json.dumps(dict(status=result["status"],rows=len(rows),rows_sha256=result["rows_sha256"],
                           counts=totals),sort_keys=True))


if __name__=="__main__":
    main()
