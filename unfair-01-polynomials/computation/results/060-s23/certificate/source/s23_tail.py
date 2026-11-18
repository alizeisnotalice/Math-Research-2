#!/usr/bin/env python3
"""Tail S={2,3}, using exact top pairs and two translated low rows."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import gcd
from pathlib import Path

from s23_entry import (
    BOUND, ONE, ZERO, ci, coeff_base, conj, extend, house, newton_base,
    plus, range_of_shift, real_sign, scale, sub,
)

TOP=((0,0),(1,0),(0,1))


def scaled_norm(P,D):
    H=[0]*17
    for i,(u,v) in enumerate(P):
        for j,(r,s) in enumerate(P):
            H[i+j]+=u*r-v*s*D
    assert all(x%4==0 for x in H)
    return [x//4 for x in H]


def disk_ascending(H):
    n=len(H)-1
    while n:
        a,b=H[n],H[0]
        if abs(a)<=abs(b):
            return n
        z=[a*H[i+1]-b*H[n-1-i] for i in range(n)]
        divisor=gcd(*z)
        H=[x//divisor for x in z]
        n-=1
    return 0


def root_code(P,D):
    for j,(p,q) in enumerate(((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))):
        value=ZERO
        for i,x in enumerate(P):
            value=plus(value,scale(p**i*q**(8-i),x))
        if real_sign(value,D)<=0 or real_sign(conj(value),D)<=0:
            return 10+j
    H=scaled_norm(P,D)
    failed=disk_ascending([h<<i for i,h in enumerate(H)])
    if failed:
        return 30+failed
    failed=disk_ascending([h<<i for i,h in enumerate(reversed(H))])
    return 50+failed if failed else 255


def select_c(r,s,g3,q3,rho4,sigma4,clo,chi,A1,B1):
    out={}
    for f in range(-9,10):
        for ea,eb in TOP:
            if B1*ea or A1*eb:
                continue
            zero_a=r+s*(eb-f)!=0
            zero_b=r+s*(ea-f)!=0
            if (zero_a and s*g3%r) or (zero_b and s*q3%r):
                continue
            ca=rho4+s*g3//r
            cb=sigma4+s*q3//r
            if zero_a and zero_b and ca!=cb:
                continue
            if zero_a:
                choices=(ca,)
            elif zero_b:
                choices=(cb,)
            else:
                choices=range(clo,chi+1)
            for c in choices:
                if clo<=c<=chi:
                    out.setdefault(c,[]).append((f,ea,eb,zero_a,zero_b))
    return out


def select_e(r,s,ea,eb,f,A7,B7,lo,hi):
    on_a=(eb-f)*(r-s*f)
    on_b=(ea-f)*(r-s*f)
    if A7!=ZERO and B7!=ZERO and on_a!=on_b:
        return ()
    if A7!=ZERO or B7!=ZERO:
        on=on_a if A7!=ZERO else on_b
        return (on//s,) if on%s==0 and lo<=on//s<=hi else ()
    return range(lo,hi+1)


def eighth_allowed(A3,B3,A8,B8,ea,eb,r,s,d,e,f):
    if A8==ZERO and B8==ZERO:
        return True
    ga1,gb1=ea-f,eb-f
    ga2,gb2=-e-f*ga1,-e-f*gb1
    for xa,xb in TOP:
        if (B3!=ZERO and xa) or (A3!=ZERO and xb):
            continue
        ga3,gb3=xa-d-e*ga1-f*ga2,xb-d-e*gb1-f*gb2
        if (A8==ZERO or r*gb2+s*gb3==0) and (B8==ZERO or r*ga2+s*ga3==0):
            return True
    return False


def one_prefix(row):
    D,L,Ta,r,Tb,s,A1,B1,g1,g2,g3,q1,q2,q3,rho4,sigma4=row
    tau=1 if D%4==1 else 2
    a,b=(Ta,r*tau),(Tb,s*tau)
    P3=[ONE,ci(L),a,b]
    Pb3=list(map(conj,P3))
    G3,Q3=(1,g1,g2,g3),(1,q1,q2,q3)
    A=(ONE,ci(A1),coeff_base(P3,G3,2),coeff_base(P3,G3,3))
    B=(ONE,ci(B1),coeff_base(Pb3,Q3,2),coeff_base(Pb3,Q3,3))
    A=(*A[:2],plus(A[2],ci(g2)),plus(A[3],ci(g3)),
       plus(coeff_base(P3+[ZERO],G3,4),ci(rho4)),ZERO)
    B=(*B[:2],plus(B[2],ci(q2)),plus(B[3],ci(q3)),
       plus(coeff_base(Pb3+[ZERO],Q3,4),ci(sigma4)),ZERO)
    s1=ci(-L)
    s2=sub(ci(L*L),scale(2,a))
    s3=plus(plus(ci(-L**3),scale(3*L,a)),scale(-3,b))
    S3=[ci(8),s1,s2,s3]
    base4=newton_base(P3,S3,4,D)
    clo,chi=range_of_shift(base4,4,D)
    cset=select_c(r,s,g3,q3,rho4,sigma4,clo,chi,A1,B1)
    counts=Counter(c_groups=len(cset))
    cases={}
    for c, top_options in sorted(cset.items()):
        g4,q4=rho4-c,sigma4-c
        G4,Q4=G3+(g4,),Q3+(q4,)
        P4,Pb4=P3+[ci(c)],Pb3+[ci(c)]
        S4=S3+[plus(base4,ci(-4*c))]
        base5=newton_base(P4,S4,5,D)
        dlo,dhi=range_of_shift(base5,5,D)
        a5,b5=coeff_base(P4+[ZERO],G4,5),coeff_base(Pb4+[ZERO],Q4,5)
        assert a5[1]==b5[1]==0 and a5[0]%2==b5[0]%2==0
        for d in range(dlo,dhi+1):
            counts["newton5_cd"]+=1
            g5,q5=-a5[0]//2-d,-b5[0]//2-d
            G5,Q5=G4+(g5,),Q4+(q5,)
            P5,Pb5=P4+[ci(d)],Pb4+[ci(d)]
            S5=S4+[plus(base5,ci(-5*d))]
            base6=newton_base(P5,S5,6,D)
            elo,ehi=range_of_shift(base6,6,D)
            if elo>ehi:continue
            P6,Pb6=P5+[ZERO],Pb5+[ZERO]
            s6_dummy=base6
            s7_dummy=newton_base(P6,S5+[s6_dummy],7,D)
            p8_dummy=newton_base(P6+[ZERO],S5+[s6_dummy,s7_dummy],8,D)
            p8_dummy=plus(p8_dummy,ci(-8))
            row6=extend(A,B,coeff_base(P6,G5,6),coeff_base(Pb6,Q5,6),6,D)
            for rho6,sigma6,A6,B6 in row6:
                admissible=[(f,ea,eb) for f,ea,eb,za,zb in top_options
                            if (not za or A6==ZERO) and (not zb or B6==ZERO)]
                if not admissible:continue
                counts["row6"]+=1
                A6s,B6s=A+(A6,),B+(B6,)
                G6,Q6=G5+(rho6,),Q5+(sigma6,)
                P7,Pb7=P6+[ZERO],Pb6+[ZERO]
                row7=extend(A6s,B6s,coeff_base(P7,G6,7),coeff_base(Pb7,Q6,7),7,D)
                for rho7,sigma7,A7,B7 in row7:
                    counts["row7"]+=1
                    A7s,B7s=A6s+(A7,),B6s+(B7,)
                    baseA8=coeff_base(P7+[ONE],G6+(rho7,),8)
                    baseB8=coeff_base(Pb7+[ONE],Q6+(sigma7,),8)
                    slopeA8=sub(ci(g2-L*g1+L*L),a)
                    slopeB8=sub(ci(q2-L*q1+L*L),conj(a))
                    for f,ea,eb in admissible:
                        for e in select_e(r,s,ea,eb,f,A7,B7,elo,ehi):
                            s7=plus(s7_dummy,ci(7*L*e-7*f))
                            if not house(s7,*BOUND[7],D):continue
                            s8=plus(plus(p8_dummy,scale(8*e,sub(a,ci(L*L)))),ci(8*L*f))
                            if not house(s8,*BOUND[8],D):continue
                            counts["newton8_top_option"]+=1
                            # f shifts the unnormalized row8 by an integer, so it does not change A8,B8 choices.
                            row8=extend(A7s,B7s,plus(baseA8,scale(e,slopeA8)),
                                         plus(baseB8,scale(e,slopeB8)),8,D)
                            for _,_,A8,B8 in row8:
                                if not eighth_allowed(A[3],B[3],A8,B8,ea,eb,r,s,d,e,f):
                                    continue
                                counts["row8_top_options"]+=1
                                mother=(D,L,Ta,r,Tb,s,c,d,e,f)
                                if mother not in cases:
                                    cases[mother]=root_code(P5+[ci(e),ci(f),ONE],D)
    return counts,cases


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("front",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    raw=a.front.read_bytes()
    rows=[tuple(json.loads(x)) for x in raw.splitlines()]
    assert len(rows)==len(set(rows))
    counts=Counter()
    cases={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for ct,res in pool.map(one_prefix,rows,chunksize=8):
            counts.update(ct)
            for key,value in res.items():
                assert key not in cases or cases[key]==value
                cases[key]=value
    summary=Counter(cases.values())
    output=b"".join(json.dumps((key,value),separators=(",",":")).encode()+b"\n"
                   for key,value in sorted(cases.items()))
    a.out.with_suffix(".rows.jsonl").write_bytes(output)
    result=dict(status="complete",front_rows=len(rows),counts=dict(counts),mothers=len(cases),
                codes=dict(sorted(summary.items())),rows_sha256=hashlib.sha256(output).hexdigest(),
                front_sha256=hashlib.sha256(raw).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
