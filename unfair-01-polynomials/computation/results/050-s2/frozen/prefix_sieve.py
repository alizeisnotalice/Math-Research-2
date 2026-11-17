#!/usr/bin/env python3
"""Exact low-prefix necessary sieve for degree-7 S={2}.

No floats, no bounded cofactor degree.  It uses the degree-7 root/Newton
necessary bounds and the known high-end consequence G3=Q3=0, then enumerates
only T,N,b and the three possible binary first rows.  Every surviving record
is a necessary low prefix, not a factorization or counterexample.
"""
from __future__ import annotations
import json, hashlib, math, sys
from fractions import Fraction as F
from pathlib import Path


def add(x,y,D): return (x[0]+y[0], x[1]+y[1])
def mul(x,y,D): return (x[0]*y[0]+x[1]*y[1]*D, x[0]*y[1]+x[1]*y[0])
def scale(x,n): return (x[0]*n,x[1]*n)
def conv(a,b,D,n):
    out=[(F(0),F(0)) for _ in range(n+1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=n: out[i+j]=add(out[i+j],mul(x,y,D),D)
    return out


def lt(x, bound, D):
    """Exact test x0+x1 sqrt(D) < bound, with x pair of Fractions."""
    p,q=x
    # Clear positive denominators: B sqrt(D) < M.
    den=p.denominator*q.denominator*bound.denominator
    A=int(p*den)
    B=int(q*den)
    C=int(bound*den)
    M=C-A
    if B==0: return 0 < M
    if B>0: return M>0 and B*B*D < M*M
    # B sqrt(D) is negative.
    return True if M>=0 else B*B*D > M*M


def unit(x,D):
    return not lt(x,F(0),D) and not lt(scale(x,-1),F(-1),D)


def abs_lt(x, bound, D):
    return lt(x,bound,D) and lt(scale(x,-1),bound,D)


def integer_value(x):
    assert x[1]==0 and x[0].denominator==1
    return x[0].numerator


def strict_integer_upper(x):
    """Largest integer n with n<x, for positive rational x."""
    return (x.numerator-1)//x.denominator


def core_tn_b():
    """Complete b-dependent core from the correct degree-7 coefficient bounds.

    b is a degree-one coefficient, so |b|<17/2.  The nonrational a is the
    degree-two coefficient, so both embeddings satisfy |a|<117/4.  The
    two s2 embeddings are b^2-T +/- sqrt(Delta), whence
    |b^2-T|+sqrt(Delta)<55/4.
    """
    ans=[]
    s2bound, abound=F(55,4),F(117,4)
    for b in range(-8,9):
        # Strict |b^2-T|<55/4, with integral T.
        for T in range(b*b-13,b*b+14):
            s2room=s2bound-abs(b*b-T)
            aroom=2*abound-abs(T)  # |T|+sqrt(Delta)<117/2.
            if s2room<=0 or aroom<=0: continue
            upper=strict_integer_upper(min(s2room*s2room,aroom*aroom))
            for D in range(1,upper+1):
                if D%4 != (T*T)%4 or math.isqrt(D)**2==D: continue
                N=(T*T-D)//4
                assert T*T-4*N==D
                ans.append((T,N,D,b))
    return ans


def build_p(T,D,b,c):
    a=(F(T,2),F(1,2)); abar=(F(T,2),F(-1,2))
    p=[(F(1),F(0)),(F(b),F(0)),a,(F(c),F(0))]
    pb=[(F(1),F(0)),(F(b),F(0)),abar,(F(c),F(0))]
    return p,pb


def main():
    out=Path(__file__).with_name('s2_low_prefix_result.json')
    tnb=core_tn_b()
    limit2,limit3=F(55,4),F(203,8)
    stage_b=[]; stage_uv=[]; survivors=[]
    for T,N,D,b in tnb:
        floor_a=(T+math.isqrt(D))//2
        x=-floor_a; y=1-T-x  # a+x and abar+y lie strictly in (0,1).
        a=(F(T,2),F(1,2)); abar=(F(T,2),F(-1,2))
        assert unit(add(a,(F(x),F(0)),D),D) and unit(add(abar,(F(y),F(0)),D),D)
        assert add(add(a,(F(x),F(0)),D),add(abar,(F(y),F(0)),D),D)==(F(1),F(0))
        s2=(F(b*b-T),F(-1))  # b^2-2a
        assert abs_lt(s2,limit2,D) and abs_lt((s2[0],-s2[1]),limit2,D)
        stage_b.append((T,N,D,b,x,y))
        for u,v in ((0,0),(1,0),(0,1)):
            g1,q1=u-b,v-b
            g2,q2=x-b*g1,y-b*q1
            # G3=Q3=0 is the high-end necessary condition for t=5.
            stage_uv.append((T,N,D,b,x,y,u,v,g1,g2,q1,q2))
            p0,pb0=build_p(T,D,b,0)
            G=[(F(1),F(0)),(F(g1),F(0)),(F(g2),F(0)),(F(0),F(0))]
            Q=[(F(1),F(0)),(F(q1),F(0)),(F(q2),F(0)),(F(0),F(0))]
            A0,B0=conv(G,p0,D,3),conv(Q,pb0,D,3)
            C0=conv(A0,B0,D,3)
            base=integer_value(C0[3])
            # Raising c by one raises C3 by two; solve C3 in {0,1} exactly.
            for eps in (0,1):
                if (eps-base)%2: continue
                c=(eps-base)//2
                if not -52<=c<=52: continue
                # s3=-b^3+3ab-3c, tested in both embeddings.
                s3=add((F(-b**3-3*c),F(0)),scale(a,3*b),D)
                if not (abs_lt(s3,limit3,D) and abs_lt((s3[0],-s3[1]),limit3,D)): continue
                p,pb=build_p(T,D,b,c)
                A,B=conv(G,p,D,3),conv(Q,pb,D,3)
                C=conv(A,B,D,3)
                if not all(unit(z,D) for z in A+B): continue
                if not all(z[1]==0 and z[0] in (0,1) for z in C): continue
                assert integer_value(C[1])==u+v and integer_value(C[2])==1 and integer_value(C[3])==eps
                survivors.append({'T':T,'N':N,'Delta':D,'b':b,'c':c,'A1':u,'B1':v,
                                  'G1':g1,'G2':g2,'G3':0,'Q1':q1,'Q2':q2,'Q3':0,
                                  'c3':eps})
    report={'status':'completed_exact_necessary_prefix_sieve','scope':'degree7 S={2}, rows 0..3 plus G3=Q3=0; not global exclusion',
            'counts':{'TN_b_complete_s2_and_a_boxes':len(tnb),'first_row_choices':len(stage_uv),'surviving_prefixes':len(survivors)},
            'core_records':[{'T':T,'N':N,'Delta':D,'b':b} for T,N,D,b in tnb],
            'survivors':survivors,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':sys.version}
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'counts':report['counts'],'out':str(out)},ensure_ascii=False))


if __name__=='__main__': main()
