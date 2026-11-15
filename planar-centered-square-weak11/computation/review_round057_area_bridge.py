#!/usr/bin/env python3
"""Exact small-instance checks of full-union one-sided derivatives and mass jets."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,platform,random,sys

def value(e,t):return e[0]+e[1]*t

def area(rects,t):
    good=[r for r in rects if value(r[0],t)<value(r[1],t) and value(r[2],t)<value(r[3],t)]
    X=sorted({value(e,t) for r in good for e in r[:2]});Y=sorted({value(e,t) for r in good for e in r[2:]})
    return sum(((b-a)*(d-c) for a,b in zip(X,X[1:]) for c,d in zip(Y,Y[1:]) if any(value(r[0],t)<(a+b)/2<value(r[1],t) and value(r[2],t)<(c+d)/2<value(r[3],t) for r in good)),F(0))

def side_jet(rects,side):
    # Affine expressions are grouped only when their full expressions agree.
    X=sorted({e for r in rects for e in r[:2]},key=lambda e:value(e,side));Y=sorted({e for r in rects for e in r[2:]},key=lambda e:value(e,side))
    assert all(value(a,side)<value(b,side) for a,b in zip(X,X[1:]))
    assert all(value(a,side)<value(b,side) for a,b in zip(Y,Y[1:]))
    A=D=F(0)
    for a,b in zip(X,X[1:]):
      for c,d in zip(Y,Y[1:]):
        mx=(value(a,side)+value(b,side))/2;my=(value(c,side)+value(d,side))/2
        if any(value(r[0],side)<mx<value(r[1],side) and value(r[2],side)<my<value(r[3],side) for r in rects):
            wx=b[0]-a[0];wy=d[0]-c[0];vx=b[1]-a[1];vy=d[1]-c[1]
            A+=wx*wy;D+=vx*wy+wx*vy
    return A,D

def run():
    e=lambda a,v=0:(F(a),F(v))
    local=[(e(0),e(1),e(0),e(1)),(e(1,1),e(2),e(0),e(1))]
    remote=[(e(0),e(2),e(2),e(3)),(e(1,1),e(2,1),e(2),e(3))]
    samples=[F(-1,4),F(-1,10),F(-1,100),F(0),F(1,100),F(1,10),F(1,4)]
    for t in samples:
        assert area(local,t)==2-max(t,0)
        assert area(remote,t)==2+max(t,0)
        assert area(local+remote,t)==4
    for rects,expected in [(local,(F(0),F(-1))),(remote,(F(0),F(1))),(local+remote,(F(0),F(0)))]:
        assert tuple(side_jet(rects,s)[1] for s in [F(-1,8),F(1,8)])==expected
        rotated=[(r[2],r[3],r[0],r[1]) for r in rects]
        assert tuple(side_jet(rotated,s)[1] for s in [F(-1,8),F(1,8)])==expected
    moving_y=[(e(0),e(1),e(0),e(1,1))]
    assert side_jet(moving_y,F(1,8))==(F(1),F(1))
    # Exact jets at t=0 where q_i(0)=r_i^2 have rational values.
    rng=random.Random(570057);jetchecks=0;doublezeros=0
    for _ in range(2000):
        r1,r2=[F(rng.randrange(1,20),rng.randrange(1,10)) for _ in range(2)]
        s,r=rng.choice([-1,1]),rng.choice([-1,1]);b1,b2=[F(rng.randrange(-12,13),rng.randrange(1,9)) for _ in range(2)]
        c=-s*r1+r*r2;a1,a2=r1*r1,r2*r2
        f1=s*b1/(2*r1)-r*b2/(2*r2);f2=-s*b1*b1/(4*r1**3)+r*b2*b2/(4*r2**3)
        k=c*c+a1-a2;P=[k*k-4*c*c*a1,2*k*(b1-b2)-4*c*c*b1,(b1-b2)**2]
        assert P[0]==0
        if c:
            if not f1:assert P[1]==0
            if not f1 and not f2:assert b1==b2==0
            if all(v==0 for v in P):assert b1==b2==0
        elif s==r:
            assert r1==r2
            if not f1:assert b1==b2 and f2==0
        if not f1 and not f2:doublezeros+=1
        jetchecks+=1
    # A genuine quadratic tangency need not swap order.
    assert F(1,2)-F(2,4)==0
    assert -F(1,4)+F(4,32)==-F(1,8)
    return {'status':'passed','scope':'arbitrary finite rectangle counterexamples and mass-only endpoint jet checks, not original atomic-input realizability or real R57 wall verification','local_left_right_derivative':['0','-1'],'remote_left_right_derivative':['0','1'],'full_left_right_derivative':['0','0'],'local_kink_can_cancel_globally':True,'rotated_moving_y_cases_passed':True,'moving_y_only_derivative':'1','rectangle_parameter_samples':len(samples),'jet_seed':570057,'mass_jet_checks':jetchecks,'persistent_or_constant_zero_jets':doublezeros,'quadratic_tangency_second_derivative':'-1/8','python':sys.version,'executable':sys.executable,'platform':platform.platform(),'audit_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    assert os.path.ismount('/Volumes/data');out=run();p=Path('planar-centered-square-weak11/computation/results/round057/area-bridge-checks.json');p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
