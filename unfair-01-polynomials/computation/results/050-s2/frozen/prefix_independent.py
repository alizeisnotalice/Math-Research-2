#!/usr/bin/env python3
"""Independent exact recheck for the corrected degree-7 S={2} prefix sieve.

This intentionally does not import scan_s2_low_prefix.py.  It enumerates
b=-8..8, T=-58..58 and Delta=1..189 directly from the two strict surd
inequalities, then uses direct coefficient formulas through degree 3.  It
compares full core and survivor *sets* against the repaired ledger.
"""
from __future__ import annotations
import hashlib, json, math, sys
from fractions import Fraction as F
from pathlib import Path


def lt_surd(p: F, q: F, D: int, L: F) -> bool:
    """p+q sqrt(D)<L by integer sign/square comparisons."""
    common=p.denominator*q.denominator*L.denominator
    A,B,C=int(p*common),int(q*common),int(L*common)
    M=C-A
    if B == 0: return M > 0
    if B > 0: return M > 0 and B*B*D < M*M
    return True if M >= 0 else B*B*D > M*M


def abs_lt(p: F, q: F, D: int, L: F) -> bool:
    return lt_surd(p,q,D,L) and lt_surd(-p,-q,D,L)


def in_unit(p: F, q: F, D: int) -> bool:
    # 0<=p+q sqrtD<=1, with equality allowed for rational rows.
    return not lt_surd(p,q,D,F(0)) and not lt_surd(-p,-q,D,F(-1))


def plus(x,y): return (x[0]+y[0],x[1]+y[1])
def times(x,y,D): return (x[0]*y[0]+x[1]*y[1]*D,x[0]*y[1]+x[1]*y[0])


def corrected_core():
    out=set()
    s2,a_bound=F(55,4),F(117,4)
    for b in range(-8,9):
        for T in range(-58,59):
            for D in range(1,190):
                if math.isqrt(D)**2==D or D%4 != T*T%4: continue
                # Two conditions are the exact complete boxes:
                # |b^2-T|+sqrtD<s2 and |T|+sqrtD<2*a_bound.
                room=s2-abs(b*b-T)
                aroom=2*a_bound-abs(T)
                if room<=0 or aroom<=0: continue
                if not (F(D)<room*room and F(D)<aroom*aroom): continue
                N=(T*T-D)//4
                out.add((T,N,D,b))
    return out


def independently_surviving(core):
    ans=set()
    L3=F(203,8)
    for T,N,D,b in core:
        floor_a=(T+math.isqrt(D))//2
        x=-floor_a
        y=1-T-x
        # a=(T+sqrtD)/2, abar=(T-sqrtD)/2.
        assert in_unit(F(T,2)+x,F(1,2),D)
        assert in_unit(F(T,2)+y,F(-1,2),D)
        for u,v in ((0,0),(1,0),(0,1)):
            g1,q1=u-b,v-b
            g2,q2=x-b*g1,y-b*q1
            for c in range(-52,53):
                # G3=Q3=0 is the independent high-end input.
                A0=(F(1),F(0)); A1=(F(u),F(0)); A2=(F(T,2)+x,F(1,2))
                B0=(F(1),F(0)); B1=(F(v),F(0)); B2=(F(T,2)+y,F(-1,2))
                A3=(F(c+b*g2,1)+F(g1*T,2),F(g1,2))
                B3=(F(c+b*q2,1)+F(q1*T,2),F(-q1,2))
                if not (in_unit(*A3,D) and in_unit(*B3,D)): continue
                s3=(F(-b**3-3*c,1)+F(3*b*T,2),F(3*b,2))
                if not (abs_lt(*s3,D,L3) and abs_lt(s3[0],-s3[1],D,L3)): continue
                # Direct degree-three convolution, rather than the other
                # program's c=0-plus-slope shortcut.
                C1=plus(times(A0,B1,D),times(A1,B0,D))
                C2=plus(plus(times(A0,B2,D),times(A1,B1,D)),times(A2,B0,D))
                C3=plus(plus(plus(times(A0,B3,D),times(A1,B2,D)),times(A2,B1,D)),times(A3,B0,D))
                if any(z[1] != 0 or z[0] not in (0,1) for z in (C1,C2,C3)): continue
                ans.add((T,N,D,b,c,u,v,g1,g2,q1,q2,int(C3[0])))
    return ans


def main():
    root=Path(__file__).parent
    repaired=root/'s2_low_prefix_result.json'
    claimed=json.loads(repaired.read_text())
    core=corrected_core()
    got_surv=independently_surviving(core)
    want_core={(z['T'],z['N'],z['Delta'],z['b']) for z in claimed['core_records']}
    want_surv={(z['T'],z['N'],z['Delta'],z['b'],z['c'],z['A1'],z['B1'],
                z['G1'],z['G2'],z['Q1'],z['Q2'],z['c3']) for z in claimed['survivors']}
    assert len(core)==4913 and core==want_core
    assert len(got_surv)==395 and got_surv==want_surv
    report={'status':'passed_independent_corrected_s2_recheck',
            'repaired_ledger_sha256':hashlib.sha256(repaired.read_bytes()).hexdigest(),
            'core_count':len(core),'survivor_count':len(got_surv),
            'core_set_equality':True,'survivor_set_equality':True,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':sys.version}
    out=root/'independent_s2_prefix_recheck.json'
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False))


if __name__=='__main__': main()
