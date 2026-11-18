#!/usr/bin/env python3
"""Exact eight-degree S={1,2} sufficient-prefix contradictions or finite entry."""

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


def candidates(D,house_two):
    tau=1 if D%4==1 else 2
    parity=1 if D%4==1 else 0
    house_double=83 if house_two else 20
    max_coef=isqrt((house_double**2-1)//(tau*tau*D))
    for n in range(-max_coef if house_two else 1,max_coef+1):
        if n==0:
            continue
        bound=house_double-isqrt((tau*n)**2*D)-1
        for T in range(-bound+((-bound-parity*n)%2),bound+1,2):
            assert (T-parity*n)%2==0
            yield (T,tau*n),n


def solve_tail(P,Pb,S,G,Q,A,B,D,r,s,counts,rows,full,state):
    cross6=product_cross(A,B,6,D)
    counts["row5"]+=1
    if not any(cross6==t for t in (ZERO,ONE)):
        return
    counts["row6_cross_pass"]+=1
    P=P[:6]+[ZERO,ZERO,ONE]
    Pb=Pb[:6]+[ZERO,ZERO,ONE]
    a6=coeff_base(P,G,6)
    b6=coeff_base(Pb,Q,6)
    assert a6[1]==b6[1]==0 and a6[0]%2==b6[0]%2==0
    g6zero,q6zero=-a6[0]//2,-b6[0]//2
    base6=newton_base(P,S,6,D)
    flo,fhi=range_of_shift(base6,6,D)
    if flo>fhi:
        return
    g_special=(r//s) if r%s==0 and -9<=r//s<=9 else None
    fixed_f_g=r*g6zero+s*G[5]
    fixed_f_q=r*q6zero+s*Q[5]
    f_fixed=fixed_f_g//r if fixed_f_g==fixed_f_q and fixed_f_g%r==0 else None
    if not full:
        return
    A6,B6=A+(ZERO,),B+(ZERO,)
    cross7=product_cross(A6,B6,7,D)
    if g_special is None:
        f_choices=[f_fixed] if f_fixed is not None and flo<=f_fixed<=fhi else []
    else:
        f_choices=range(flo,fhi+1)
    for f in f_choices:
        P[6],Pb[6]=ci(f),ci(f)
        g6,q6=g6zero-f,q6zero-f
        S6=[*S,plus(base6,ci(-6*f))]
        base7=newton_base(P,S6,7,D)
        gl,gh=range_of_shift(base7,7,D)
        if f!=f_fixed or cross7 not in (ZERO,ONE):
            g_choices=[g_special] if g_special is not None and gl<=g_special<=gh else []
        else:
            g_choices=range(gl,gh+1)
        for g in g_choices:
            P[7],Pb[7]=ci(g),ci(g)
            S7=[*S6,plus(base7,ci(-7*g))]
            s8=plus(newton_base(P,S7,8,D),ci(-8))
            if not house(s8,*BOUND[8],D):
                continue
            G6,Q6=G+(g6,),Q+(q6,)
            A7base,B7base=coeff_base(P,G6,7),coeff_base(Pb,Q6,7)
            if g==g_special:
                row7=extend(A6,B6,A7base,B7base,7,D)
            else:
                assert A7base[1]==B7base[1]==0 and A7base[0]%2==B7base[0]%2==0
                row7=[(-A7base[0]//2,-B7base[0]//2,ZERO,ZERO)]
            for g7,q7,ak,bk in row7:
                counts["row7"]+=1
                A7,B7=A6+(ak,),B6+(bk,)
                u8,v8=coeff_base(P,G6+(g7,),8),coeff_base(Pb,Q6+(q7,),8)
                accepted=list(extend(A7,B7,u8,v8,8,D))
                if accepted:
                    counts["row8_states"]+=len(accepted)
                    row=(*state,f,g)
                    if row not in rows:
                        rows[row]=root_code(P,D)


def one_domain(job):
    D,full=job
    counts=Counter()
    rows={}
    initial_b=list(candidates(D,True))
    if not initial_b:
        return counts,rows
    for a,r in candidates(D,False):
        a_bar=conj(a)
        g1=-floor_quad(a,D)
        q1=1-a[0]-g1
        A1=plus(a,ci(g1))
        B1=plus(a_bar,ci(q1))
        assert plus(A1,B1)==ONE
        A=(ONE,A1)
        B=(ONE,B1)
        aa=mult(a,a,D)
        for b,s in initial_b:
            counts["raw_core"]+=1
            s2=sub(aa,scale(2,b))
            if not house(s2,*BOUND[2],D):
                continue
            counts["newton2_core"]+=1
            b_bar=conj(b)
            row2=list(extend(A,B,plus(b,scale(g1,a)),plus(b_bar,scale(q1,a_bar)),2,D,True))
            assert len(row2)<=1
            if not row2:
                continue
            g2,q2,A2,B2=row2[0]
            A_two,B_two=A+(A2,),B+(B2,)
            counts["row2"]+=1
            row3=list(extend(A_two,B_two,plus(scale(g1,b),scale(g2,a)),
                             plus(scale(q1,b_bar),scale(q2,a_bar)),3,D,True))
            assert len(row3)<=1
            if not row3:
                continue
            counts["row3_normalized"]+=1
            rho3,sigma3,A3,B3=row3[0]
            A_three,B_three=A_two+(A3,),B_two+(B3,)
            P=[ONE,a,b]
            S=[ci(8),neg(a),s2]
            base3=newton_base(P,S,3,D)
            clo,chi=range_of_shift(base3,3,D)
            d_slope=sub(ci(g1-q1),sub(a,a_bar))
            assert real_sign(plus(d_slope,ci(2)),D)<0
            for c in range(clo,chi+1):
                counts["newton3"]+=1
                G_three,Q_three=(1,g1,g2,rho3-c),(1,q1,q2,sigma3-c)
                P3=P+[ci(c)]
                Pb3=list(map(conj,P3))
                S3=[*S,plus(base3,ci(-3*c))]
                x4=plus(plus(scale(g1,ci(c)),scale(g2,b)),scale(G_three[3],a))
                y4=plus(plus(scale(q1,ci(c)),scale(q2,b_bar)),scale(Q_three[3],a_bar))
                row4=list(extend(A_three,B_three,x4,y4,4,D,True))
                assert len(row4)<=1
                if not row4:
                    continue
                counts["row4_normalized"]+=1
                rho4,sigma4,A4,B4=row4[0]
                if s*(rho4-sigma4)%r!=0:
                    continue
                counts["row4_difference_divisibility"]+=1
                A_four,B_four=A_three+(A4,),B_three+(B4,)
                base4=newton_base(P3,S3,4,D)
                dlo,dhi=range_of_shift(base4,4,D)
                if dlo>dhi:
                    continue
                difference=plus(plus(ci(c*(g2-q2)-s*(rho4-sigma4)//r),
                                     sub(scale(G_three[3],b),scale(Q_three[3],b_bar))),
                                sub(scale(rho4,a),scale(sigma4,a_bar)))
                d=min_d_with_bound(difference,d_slope,1,dlo,dhi,D)
                if d>dhi or real_sign(plus(plus(difference,scale(d,d_slope)),ci(1)),D)<0:
                    continue
                counts["row5_difference"]+=1
                g4,q4=rho4-d,sigma4-d
                if s*g4%r or s*q4%r:
                    continue
                g5,q5=-s*g4//r,-s*q4//r
                G,Q=G_three+(g4,g5),Q_three+(q4,q5)
                P4,Pb4=P3+[ci(d)],Pb3+[ci(d)]
                x5,y5=coeff_base(P4+[ZERO],G,5),coeff_base(Pb4+[ZERO],Q,5)
                x5,y5=plus(x5,ci(g5)),plus(y5,ci(q5))
                c5=plus(plus(x5,y5),product_cross(A_four,B_four,5,D))
                assert c5[1]==0 and c5[0]%2==0
                e=-(c5[0]//4)
                if not -92<=e<=92:
                    continue
                A5,B5=plus(x5,ci(e)),plus(y5,ci(e))
                if not unit(A5,D) or not unit(B5,D):
                    continue
                assert plus(c5,ci(2*e)) in (ZERO,ONE)
                S4=[*S3,plus(base4,ci(-4*d))]
                base5=newton_base(P4,S4,5,D)
                s5=plus(base5,ci(-5*e))
                if not house(s5,*BOUND[5],D):
                    continue
                state=(D,a[0],r,b[0],s,c,d,e)
                solve_tail(P4+[ci(e)],Pb4+[ci(e)],[*S4,s5],G,Q,A_four+(A5,),B_four+(B5,),
                           D,r,s,counts,rows,full,state)
    counts["mothers"]=len(rows)
    counts.update(rows.values())
    return counts,rows


def degree8_supports():
    banned={(1,2),(1,3),(1,4),(2,1),(3,1),(4,1),(2,4),(4,2)}
    def rev(s):
        return tuple(8-i for i in s[::-1])
    keep=[]
    for mask in range(1,128):
        s=tuple(i+1 for i in range(7) if mask&(1<<i))
        r,t=s[0],8-s[-1]
        if r!=t and (r,t) not in banned:
            keep.append(s)
    return sorted({min(s,rev(s)) for s in keep},key=lambda s:(len(s),s))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--full",action="store_true")
    args=ap.parse_args()
    assert not args.out.exists()
    domains=[D for D in range(2,400) if all(D%(p*p)!=0 for p in range(2,isqrt(D)+1))
             and any(candidates(D,False))]
    all_rows={}
    counts=Counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for local,rows in pool.map(one_domain,[(D,args.full) for D in domains],chunksize=1):
            counts.update(local)
            assert all_rows.keys().isdisjoint(rows)
            all_rows.update(rows)
    cases=sorted((row,code) for row,code in all_rows.items())
    payload=b"".join(json.dumps(rec,separators=(",",":")).encode()+b"\n" for rec in cases)
    args.out.with_suffix(".rows.jsonl").write_bytes(payload)
    summary=dict(status="complete",full=args.full,domains=len(domains),counts=dict(counts),
                 cases_sha256=hashlib.sha256(payload).hexdigest(),reversal_orbits=degree8_supports(),
                 initial_domain_range=[min(domains),max(domains)])
    args.out.write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print(json.dumps(summary,sort_keys=True))


if __name__=="__main__":
    main()
