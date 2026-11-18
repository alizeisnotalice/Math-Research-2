#!/usr/bin/env python3
"""Exact eight-degree S={2,5} sufficient-prefix contradictions or finite entry."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import gcd, isqrt
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


def domain_a(D):
    t=int(D%4==1)
    tau=1 if t else 2
    for r in range(1,isqrt(288//(tau*tau*D))+1):
        R=16-isqrt((tau*r)**2*D)
        for T in range(-R+(t*r+R)%2,R+1,2):
            yield (T,tau*r),r


def possible_b(a,r,c,D):
    t=int(D%4==1)
    tau=1 if t else 2
    for s in range(-r,r+1):
        if not s or r%s:
            continue
        w=(t*s,tau*s)
        if s==r*c:
            for k in range(-25,26):
                b=plus(scale(c,a),ci(k))
                assert b[1]==tau*s
                yield b,s,True
        else:
            for u in sorted({r//s-floor_quad(w,D),r//s-floor_quad(conj(w),D)}):
                yield plus(w,ci(u)),s,False


def top_choices(ea,eb,f,a4,b4,exception):
    if a4!=ZERO and eb!=f:
        return False
    if b4!=ZERO and ea!=f:
        return False
    return exception or f==0


def tail_e(r,s,rho,sigma,ea,eb,f,exception):
    a1,b1=ea-f,eb-f
    # R3(e)=floor-shift(b) - f^2(f-eps) + e(2f-eps).
    const_a=r+s*(rho+f*f*a1)
    const_b=r+s*(sigma+f*f*b1)
    slope_a=s*(f-a1)
    slope_b=s*(f-b1)
    if not exception:
        return (0,) if const_a==0 or const_b==0 else ()
    if (slope_a==0 and const_a==0) or (slope_b==0 and const_b==0):
        return range(-41,42)
    result=set()
    for const,slope in ((const_a,slope_a),(const_b,slope_b)):
        if slope and const%slope==0:
            k=-const//slope
            if -41<=k<=41:
                result.add(k)
    return sorted(result)


def required_zero_choices(base,D,need_zero):
    if not need_zero:
        return shifts(base,D)
    if base[1]:
        return ()
    assert base[0]%2==0
    return ((-base[0]//2,ZERO),)


def reverse_coupled(P,D,G,Q,A,B,low_a,low_b):
    R=P[::-1]
    S=list(map(conj,R))
    states=[(G,Q,A,B)]
    for k in range(4,9):
        triggers=(low_b[k]!=ZERO,low_a[k]!=ZERO)
        next_layer=[]
        for U,V,C,E in states:
            cx=product_cross(C,E,k,D)
            for ui,ck in required_zero_choices(coeff_base(R,U,k),D,triggers[0]):
                for vi,ek in required_zero_choices(coeff_base(S,V,k),D,triggers[1]):
                    if plus(cx,plus(ck,ek)) in (ZERO,ONE):
                        next_layer.append((U+(ui,),V+(vi,),C+(ck,),E+(ek,)))
        if not next_layer:
            return False
        states=next_layer
    return True


def one_domain(D):
    cases={}
    counts=Counter()
    for a,r in domain_a(D):
        counts["a_cores"]+=1
        assert house(a,17,2,D)
        ab=conj(a)
        g2,q2=-floor_quad(a,D),-floor_quad(ab,D)
        a2,b2=plus(a,ci(g2)),plus(ab,ci(q2))
        assert plus(a2,b2)==ONE and (g2 or q2)
        A,B=(ONE,ZERO,a2,ZERO),(ONE,ZERO,b2,ZERO)
        row4=list(extend(A,B,scale(g2,a),scale(q2,ab),4,D,True))
        assert len(row4)<=1
        for rho4,sigma4,a4,b4 in row4:
            counts["a_row4"]+=1
            assert (a4!=ZERO)==bool(g2) and (b4!=ZERO)==bool(q2)
            A4,B4=A+(a4,),B+(b4,)
            for c in range(-10,11):
                G3,Q3=(1,0,g2,-c),(1,0,q2,-c)
                S3=[ci(8),ZERO,scale(-2,a),ci(-3*c)]
                ns4=scale(2,mult(a,a,D))
                dlo,dhi=range_of_shift(ns4,4,D)
                for b,s,exception in possible_b(a,r,c,D):
                    counts["initial_b"]+=1
                    bb=conj(b)
                    if not house(b,185,2,D):
                        continue
                    s5=scale(5,sub(scale(c,a),b))
                    if not house(s5,*BOUND[5],D):
                        continue
                    counts["newton5_ab"]+=1
                    rho_b,sigma_b=-floor_quad(b,D),-floor_quad(bb,D)
                    ar3,br3=plus(b,ci(rho_b)),plus(bb,ci(sigma_b))
                    assert plus(ar3,br3)==ONE
                    a5start=plus(b,plus(ci(c*g2),scale(-c,a)))
                    b5start=plus(bb,plus(ci(c*q2),scale(-c,ab)))
                    row5=list(extend(A4,B4,a5start,b5start,5,D))
                    for ea,eb in TOP:
                        for f in range(2):
                            if not top_choices(ea,eb,f,a4,b4,exception):
                                continue
                            a1,b1=ea-f,eb-f
                            for e in tail_e(r,s,rho_b,sigma_b,ea,eb,f,exception):
                                counts["top_options"]+=1
                                u2,v2=-e-f*a1,-e-f*b1
                                u3,v3=rho_b-e*a1-f*u2,sigma_b-e*b1-f*v2
                                za,zb=r+s*v3!=0,r+s*u3!=0
                                assert not (za and zb)
                                if za:
                                    ds=(rho4,)
                                elif zb:
                                    ds=(sigma4,)
                                else:
                                    ds=range(dlo,dhi+1)
                                Gtop,Qtop=(1,a1,u2,u3),(1,b1,v2,v3)
                                Atop,Btop=(ONE,ci(ea),ZERO,ar3),(ONE,ci(eb),ZERO,br3)
                                for d in ds:
                                    if not dlo<=d<=dhi:
                                        continue
                                    counts["newton4_top_d"]+=1
                                    P=[ONE,ZERO,a,ci(c),ci(d),b,ci(e),ci(f),ONE]
                                    Pb=list(map(conj,P))
                                    S5=[*S3,plus(ns4,ci(-4*d)),s5]
                                    s6=plus(newton_base(P,S5,6,D),ci(-6*e))
                                    if not house(s6,*BOUND[6],D):
                                        continue
                                    s7=plus(newton_base(P,S5+[s6],7,D),ci(-7*f))
                                    if not house(s7,*BOUND[7],D):
                                        continue
                                    s8=plus(newton_base(P,S5+[s6,s7],8,D),ci(-8))
                                    if not house(s8,*BOUND[8],D):
                                        continue
                                    counts["newton8_top_d"]+=1
                                    G4,Q4=G3+(rho4-d,),Q3+(sigma4-d,)
                                    for g5,q5,a5,b5 in row5:
                                        G5,Q5=G4+(g5,),Q4+(q5,)
                                        A5,B5=A4+(a5,),B4+(b5,)
                                        for g6,q6,a6,b6 in extend(A5,B5,coeff_base(P,G5,6),
                                                                    coeff_base(Pb,Q5,6),6,D,True):
                                            if (za and a6!=ZERO) or (zb and b6!=ZERO):
                                                continue
                                            counts["row6"]+=1
                                            G6,Q6=G5+(g6,),Q5+(q6,)
                                            A6,B6=A5+(a6,),B5+(b6,)
                                            for g7,q7,a7,b7 in extend(A6,B6,coeff_base(P,G6,7),
                                                                         coeff_base(Pb,Q6,7),7,D):
                                                counts["row7"]+=1
                                                G7,Q7=G6+(g7,),Q6+(q7,)
                                                A7,B7=A6+(a7,),B6+(b7,)
                                                for _,_,a8,b8 in extend(A7,B7,coeff_base(P,G7,8),
                                                                          coeff_base(Pb,Q7,8),8,D):
                                                    counts["row8"]+=1
                                                    if not reverse_coupled(P,D,Gtop,Qtop,Atop,Btop,A7+(a8,),B7+(b8,)):
                                                        continue
                                                    counts["coupled_row8"]+=1
                                                    mother=(D,a[0],r,c,d,b[0],s,e,f)
                                                    if mother not in cases:
                                                        cases[mother]=root_code(P,D)
                                                        counts["mothers_exception" if exception else "mothers_generic"]+=1
    counts["mothers"]=len(cases)
    counts.update(cases.values())
    return counts,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    args=ap.parse_args()
    assert not args.out.exists()
    domains=[D for D in range(2,289) if all(D%(p*p)!=0 for p in range(2,isqrt(D)+1))
             and any(domain_a(D))]
    totals=Counter()
    all_rows={}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for counts,rows in pool.map(one_domain,domains,chunksize=1):
            totals.update(counts)
            assert all_rows.keys().isdisjoint(rows.keys())
            all_rows.update(rows)
    lines=sorted(all_rows.items())
    payload=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in lines)
    args.out.with_suffix(".rows.jsonl").write_bytes(payload)
    summary=dict(status="complete",domains=len(domains),counts=dict(totals),
                 rows_sha256=hashlib.sha256(payload).hexdigest())
    args.out.write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print(json.dumps(summary,sort_keys=True))


if __name__=="__main__":
    main()
