#!/usr/bin/env python3
"""Exact eight-degree S={2,3,4,5} sufficient-prefix contradictions or finite entry."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import gcd, isqrt
from functools import lru_cache
from pathlib import Path

ZERO, ONE = (0,0), (2,0)
BOUND = {1:(10,1), 2:(17,1), 3:(65,2), 4:(257,4), 5:(1025,8),
         6:(4097,16), 7:(16385,32), 8:(65537,64)}


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


def screen_step(a):
    degree=len(a)-1
    while degree:
        leading,constant=a[0],a[degree]
        if abs(leading)<=abs(constant):
            return False
        out=[leading*a[i]-constant*a[degree-i] for i in range(degree)]
        divisor=gcd(*out)
        a=[x//divisor for x in out]
        degree-=1
    return True


def root_code(P,D):
    for p,q in ((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6)):
        value=ZERO
        for j,x in enumerate(P):
            value=plus(value,scale(p**j*q**(8-j),x))
        if real_sign(value,D)<=0 or real_sign(conj(value),D)<=0:
            return "root_positive_reject"
    H=[0]*17
    for i,(u,v) in enumerate(P):
        for j,(r,s) in enumerate(P):
            H[i+j]+=u*r-v*s*D
    assert all(h%4==0 for h in H)
    H=[h//4 for h in H]
    if not screen_step([H[j]<<j for j in range(16,-1,-1)]):
        return "root_outer_reject"
    if not screen_step([H[j]<<(16-j) for j in range(17)]):
        return "root_inner_reject"
    return "root_pass"

TOP=((0,0),(1,0),(0,1))


def forced_choices(x,D,force):
    if not force:return shifts(x,D)
    if x[1]:return ()
    assert x[0]%2==0
    return ((-x[0]//2,ZERO),)


def reverse_coupled(P,D,U,V,A,B,low_a,low_b):
    R=P[::-1]
    S=list(map(conj,R))
    states=[(U,V,A,B)]
    for k in range(4,9):
        nxt=[]
        for U,V,A,B in states:
            cross=product_cross(A,B,k,D)
            for ug,ak in forced_choices(coeff_base(R,U,k),D,low_b[k]!=ZERO):
                for vg,bk in forced_choices(coeff_base(S,V,k),D,low_a[k]!=ZERO):
                    if plus(cross,plus(ak,bk)) in (ZERO,ONE):
                        nxt.append((U+(ug,),V+(vg,),A+(ak,),B+(bk,)))
        if not nxt:return False
        states=nxt
    return True



def squarefree(D):
    return all(D % (p*p) for p in range(2, isqrt(D)+1))


def a_cores(D,L):
    eps=int(D%4==1)
    tau=1 if eps else 2
    for r in range(1,isqrt(288//(tau*tau*D))+1):
        radius=16-isqrt((tau*r)**2*D)
        centre=L*L
        lo,hi=centre-radius,centre+radius
        for trace in range(lo+(eps*r-lo)%2,hi+1,2):
            a=(trace,tau*r)
            if house(a,83,2,D):
                yield a,r


def shift_interval(base,delta,num,den,D):
    """All z with both embeddings of base+delta*z strictly bounded."""
    if delta<0:
        base,delta=neg(base),-delta
    if not delta:
        return (-120,120) if house(base,num,den,D) else (1,0)
    u,v=base
    R=(2*num-isqrt((den*v)**2*D)-1)//den
    return -((R+u)//(2*delta)),(R-u)//(2*delta)


def tail_e(f):
    return range(max(-41,(f*f-16+1)//2),min(41,(f*f+16)//2)+1)


@lru_cache
def tail_boxes(D,L,r):
    """Group exact necessary (t,j,e,f) choices by the integral b-La."""
    tau=1 if D%4==1 else 2
    t_max=isqrt(4224//(9*tau*tau*D))
    Bmin,Bmax=-((32+L**3)//3),(32-L**3)//3
    boxes={}
    for t in range(-t_max,t_max+1):
        if not t:continue
        for j in range(-10,10):
            if not j:continue
            q,rem=divmod(t*(1-L*j),r)
            if r%t==0 and -9<=j<=9:
                e=L*(r//t)
                if e in tail_e(j):
                    for B in range(Bmin,Bmax+1):
                        if rem or B!=q:
                            boxes.setdefault(B,{}).setdefault((t,j),[]).append((e,j,False))
            if not rem and Bmin<=q<=Bmax:
                for f in (j,j+1):
                    if -9<=f<=9:
                        target=boxes.setdefault(q,{}).setdefault((t,j),[])
                        target.extend((e,f,True) for e in tail_e(f))
    return boxes


def d_options(D,L,r,t,j,e,f,x,y,dbase):
    U1,V1=x-f,y-f
    U2,V2=-e-f*U1,-e-f*V1
    translations=set()
    for k1,k2,d0 in ((U1,U2,dbase),(V1,V2,conj(dbase))):
        h,remainder=divmod(r*(1+L*k1),t)
        if not remainder:
            translations.add(h+(j-f)*k2-e*k1-floor_quad(d0,D))
    for z in sorted(translations):
        d=plus(dbase,ci(z))
        S3=plus(ci(-f**3+3*f*e),scale(-3,d))
        if not house(S3,*BOUND[3],D):continue
        # |d|<185/2 follows independently from the coefficient estimate.
        if not house(d,185,2,D):continue
        db=conj(d)
        p,q=-floor_quad(d,D),-floor_quad(db,D)
        U3,V3=p-e*U1-f*U2,q-e*V1-f*V2
        highA=r*(1+L*U1)+t*(j*U2+U3)
        highB=r*(1+L*V1)+t*(j*V2+V3)
        assert not (highA and highB)
        X=(ONE,ci(x),ZERO,plus(d,ci(p)))
        Y=(ONE,ci(y),ZERO,plus(db,ci(q)))
        assert plus(X[3],Y[3])==ONE
        yield d,db,(1,U1,U2,U3),(1,V1,V2,V3),highA,highB,X,Y


def one_domain(D):
    totals=Counter()
    cases={}
    eps=int(D%4==1)
    tau=1 if eps else 2
    for L in range(-9,10):
        if not L:continue
        for a,r in a_cores(D,L):
            totals["a_cores"]+=1
            ab=conj(a)
            p,q=floor_quad(a,D),floor_quad(ab,D)
            g2,q2=L*L-p,L*L-q
            A2,B2=plus(a,ci(-p)),plus(ab,ci(-q))
            assert plus(A2,B2)==ONE
            A3,B3=(ONE,ZERO,A2,ZERO),(ONE,ZERO,B2,ZERO)
            ns2=sub(ci(L*L),scale(2,a))
            a4base=plus(ci(L**4),scale(2,mult(a,a,D)))
            for B,geometry in tail_boxes(D,L,r).items():
                b=plus(scale(L,a),ci(B))
                if not house(b,185,2,D):continue
                bb=conj(b)
                g3,q3=-B-L*g2,-B-L*q2
                G3,Q3=(1,-L,g2,g3),(1,-L,q2,q3)
                P3=(ONE,ci(L),a,b)
                Pb3=(ONE,ci(L),ab,bb)
                S3=(ci(8),ci(-L),ns2,ci(-L**3-3*B))
                a5base=plus(ci(-L**5-5*L*L*B),scale(5*B,a))
                for (t,j),tails in geometry.items():
                    s=t*j
                    c0=(eps*s,tau*s)
                    cb0=conj(c0)
                    ns4=sub(plus(a4base,ci(4*L*B)),scale(4,c0))
                    clo,chi=shift_interval(ns4,-4,*BOUND[4],D)
                    cbox=shift_interval(c0,1,1921,16,D)
                    clo,chi=max(clo,cbox[0]),min(chi,cbox[1])
                    if clo>chi:continue
                    totals["core_LaBtj"]+=1
                    for rho,sigma,a4,b4 in extend(A3,B3,coeff_base(P3+(c0,),G3,4),
                                                   coeff_base(Pb3+(cb0,),Q3,4),4,D,True):
                        totals["row4_normalized"]+=1
                        assert (a4!=ZERO)==bool(s-r*p) and (b4!=ZERO)==bool(s-r*q)
                        A4,B4=A3+(a4,),B3+(b4,)
                        for e,f,exception in tails:
                            dbase=(eps*t,tau*t)
                            for x,y in TOP:
                                if (a4!=ZERO and y!=f-j) or (b4!=ZERO and x!=f-j):continue
                                for d,db,U,V,highA,highB,X,Y in d_options(D,L,r,t,j,e,f,x,y,dbase):
                                    totals["high_options"]+=1
                                    lb,hb=clo,chi
                                    if highB:
                                        n=s*g2-L*t
                                        if n%r:continue
                                        lb=hb=rho+L*g3+n//r
                                    if highA:
                                        n=s*q2-L*t
                                        if n%r:continue
                                        lb=hb=sigma+L*q3+n//r
                                    lb,hb=max(lb,clo),min(hb,chi)
                                    ns5=sub(plus(a5base,scale(5*L,c0)),scale(5,d))
                                    l5,h5=shift_interval(ns5,5*L,*BOUND[5],D)
                                    lb,hb=max(lb,l5),min(hb,h5)
                                    rs4=sub(plus(ci(f**4-4*f*f*e+2*e*e),scale(4*f,d)),scale(4,c0))
                                    rl4,rh4=shift_interval(rs4,-4,*BOUND[4],D)
                                    lb,hb=max(lb,rl4),min(hb,rh4)
                                    rs5=sub(plus(plus(ci(-f**5+5*f**3*e-5*f*e*e),
                                                        scale(5*(e-f*f),d)),scale(5*f,c0)),scale(5,b))
                                    rl5,rh5=shift_interval(rs5,5*f,*BOUND[5],D)
                                    lb,hb=max(lb,rl5),min(hb,rh5)
                                    for z in range(lb,hb+1):
                                        c,cb=plus(c0,ci(z)),plus(cb0,ci(z))
                                        P=P3+(c,d,ci(e),ci(f),ONE)
                                        Pb=Pb3+(cb,db,ci(e),ci(f),ONE)
                                        S=[*S3,plus(ns4,ci(-4*z)),plus(ns5,ci(5*L*z))]
                                        for k in range(6,9):
                                            sk=plus(newton_base(P,S,k,D),scale(-k,P[k]))
                                            if not house(sk,*BOUND[k],D):break
                                            S.append(sk)
                                        if len(S)!=9:continue
                                        R=P[::-1]
                                        SR=[ci(8),ci(-f),ci(f*f-2*e),
                                            sub(ci(-f**3+3*f*e),scale(3,d)),
                                            plus(rs4,ci(-4*z)),plus(rs5,ci(5*f*z))]
                                        for k in range(6,9):
                                            sk=plus(newton_base(R,SR,k,D),scale(-k,R[k]))
                                            if not house(sk,*BOUND[k],D):break
                                            SR.append(sk)
                                        if len(SR)!=9:continue
                                        totals["newton8_both"]+=1
                                        G4,Q4=G3+(rho-z,),Q3+(sigma-z,)
                                        states=[(G4,Q4,A4,B4)]
                                        for k in range(5,9):
                                            nxt=[]
                                            for G,Q,A,Bp in states:
                                                for g,h,ak,bk in extend(A,Bp,coeff_base(P,G,k),coeff_base(Pb,Q,k),k,D):
                                                    if k==6 and ((highB and ak!=ZERO) or (highA and bk!=ZERO)):
                                                        continue
                                                    nxt.append((G+(g,),Q+(h,),A+(ak,),Bp+(bk,)))
                                            if not nxt:
                                                states=[]
                                                break
                                            states=nxt
                                            totals[f"row{k}"]+=len(nxt)
                                        for _,_,A,Bp in states:
                                            if not reverse_coupled(P,D,U,V,X,Y,A,Bp):continue
                                            totals["coupled_row8"]+=1
                                            row=(D,L,a[0],r,B,c[0],s,d[0],t,e,f)
                                            if row not in cases:
                                                cases[row]=root_code(P,D)
                                                totals["mothers_exception" if exception else "mothers_generic"]+=1
    totals["mothers"]=len(cases)
    totals.update(cases.values())
    return totals,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--domain",type=int,action="append")
    a=ap.parse_args()
    assert not a.out.exists() and not a.out.with_suffix(".rows.jsonl").exists()
    full=[D for D in range(2,289) if squarefree(D) and any(any(a_cores(D,L)) for L in (1,2))]
    domains=full if a.domain is None else sorted(set(a.domain))
    assert set(domains)<=set(full)
    counts=Counter()
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for ct,local in pool.map(one_domain,domains,chunksize=1):
            counts.update(ct)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
    payload=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in sorted(rows.items()))
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    result=dict(status="complete" if domains==full else "partial_domains",domains=len(domains),
                domain_list=domains,counts=dict(counts),rows_sha256=hashlib.sha256(payload).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
