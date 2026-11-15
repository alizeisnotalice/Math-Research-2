#!/usr/bin/env python3
"""Directed integer whole-cell source verification. No MATLAB or LP launch.

All interval operations round outwards on 2**BITS; inputs are exact JSON
numbers (decimal spelling) or rational strings. A PASS covers EVERY D4
representative cell and hence the whole support. No D-minimum claim is made.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse, hashlib, json, os, sys, time

BITS=160
S=1<<BITS

def ceildiv(a,b):
    assert b>0
    return -((-a)//b)
def iv(x):
    x=F(x)
    return (x.numerator*S//x.denominator,ceildiv(x.numerator*S,x.denominator))
def plus(a,b): return (a[0]+b[0],a[1]+b[1])
def minus(a,b): return (a[0]-b[1],a[1]-b[0])
def times(a,b):
    q=[a[i]*b[j] for i in (0,1) for j in (0,1)]
    return (min(q)//S,ceildiv(max(q),S))
def scale(a,q):
    q=F(q)
    if q<0:return scale((-a[1],-a[0]),-q)
    return (a[0]*q.numerator//q.denominator,ceildiv(a[1]*q.numerator,q.denominator))
def divide(a,b):
    assert b[0]>0
    reciprocal=(S*S//b[1],ceildiv(S*S,b[0]))
    return times(a,reciprocal)
def square(a):
    q=[a[0]*a[0],a[1]*a[1]]
    lo=0 if a[0]<=0<=a[1] else min(q)
    return (lo//S,ceildiv(max(q),S))
def power(a,n):
    assert n>=0
    out=(S,S)
    for _ in range(n):out=times(out,a)
    return out
ZERO=(0,0)
ONE=(S,S)

@lru_cache(maxsize=None)
def atan_small(x):
    assert 0<=x<=F(1,2)
    z=iv(x);sq=square(z);p=z;out=ZERO
    # 88 terms: with x<=1/2, tail <2^-176/177, below fixed grid.
    for n in range(88):
        term=scale(p,F(1,2*n+1))
        out=plus(out,term) if n%2==0 else minus(out,term)
        p=times(p,sq)
    rem=ceildiv(p[1],177)
    return (out[0]-rem,out[1]+rem)

PI=minus(scale(atan_small(F(1,5)),16),scale(atan_small(F(1,239)),4))
assert PI[0]>3*S and PI[1]*7<22*S
# Machin identity tangent and quadrant checked with exact rational arithmetic.
t2=2*F(1,5)/(1-F(1,5)**2);t4=2*t2/(1-t2*t2)
assert (t4-F(1,239))/(1+t4*F(1,239))==1

@lru_cache(maxsize=None)
def atan_pos(x):
    x=F(x);assert x>=0
    if x>1:return minus(scale(PI,F(1,2)),atan_pos(1/x))
    if x>F(1,2):return minus(scale(PI,F(1,4)),atan_small((1-x)/(1+x)))
    return atan_small(x)

def atan_signed(x):
    if x>=0:return atan_pos(x)
    a=atan_pos(-x);return (-a[1],-a[0])

@lru_cache(maxsize=None)
def boxprob_point(distance,delta,h):
    # Difference of two arctans has no quadrant ambiguity; either sign allowed.
    a=atan_signed((delta/2-distance)/h)
    b=atan_signed((-delta/2-distance)/h)
    out=divide(minus(a,b),PI)
    return (max(0,out[0]),min(S,out[1]))

@lru_cache(maxsize=None)
def cauchyjet(lo,hi,h):
    x=(iv(lo)[0],iv(hi)[1]);sq=square(x);hsq=iv(h*h)
    den=plus(hsq,sq)
    k=divide(iv(h),times(PI,den))
    kpp=divide(scale(minus(scale(sq,3),hsq),2*h),times(PI,power(den,3)))
    return k,kpp

@lru_cache(maxsize=None)
def kprime(lo,hi,h):
    x=(iv(lo)[0],iv(hi)[1]);den=plus(iv(h*h),square(x))
    return divide(scale(x,-2*h),times(PI,square(den)))

@lru_cache(maxsize=None)
def boxjet(lo,hi,center,delta,h):
    l,r=lo-center,hi-center
    dmax=max(abs(l),abs(r));dmin=F(0) if l<=0<=r else min(abs(l),abs(r))
    small=boxprob_point(dmax,delta,h);large=boxprob_point(dmin,delta,h)
    prob=(small[0],large[1])
    second=minus(kprime(l+delta/2,r+delta/2,h),kprime(l-delta/2,r-delta/2,h))
    return prob,second

def orbits(N,delta):
    out=[]
    for i in range(N+1):
        for j in range(i+1):
            x,y=i*delta,j*delta
            pts={(sx*x,sy*y) for sx in (-1,1) for sy in (-1,1)}
            pts|={(b,a) for a,b in pts}
            out.append(((x,y),sorted(pts)))
    return out

class Model:
    def __init__(self,data,repair):
        self.h=F(data['h']);self.delta=F(data['delta']);self.N=int(data['N'])
        self.lam=F(data['source_origin_mass'])+repair
        self.a=list(map(F,data['density_coefficients']))
        self.oo=orbits(self.N,self.delta)
        assert self.h>0 and self.delta>0 and self.N>=1
        assert len(self.a)==len(self.oo) and self.lam>=0 and min(self.a)>=0
        if 'orbit_representatives' in data:
            assert [tuple(map(F,p)) for p in data['orbit_representatives']]==[r for r,_ in self.oo]
        if 'orbit_sizes' in data:assert list(map(int,data['orbit_sizes']))==[len(p) for _,p in self.oo]
        # Coefficients belong to physical cells. Distinct orbit supports are disjoint.
        self.physical=[(x,y,a) for a,(_,pts) in zip(self.a,self.oo) if a for x,y in pts]
        assert sum(len(p) for _,p in self.oo)==(2*self.N+1)**2
        self.mass=self.lam+self.delta*self.delta*sum(a*len(p) for a,(_,p) in zip(self.a,self.oo))
        self.point_cache={}
    def gpoint(self,x,y):
        key=(x,y)
        if key in self.point_cache:return self.point_cache[key]
        kx=cauchyjet(x,x,self.h)[0];ky=cauchyjet(y,y,self.h)[0]
        value=scale(times(kx,ky),self.lam)
        for u,v,a in self.physical:
            value=plus(value,scale(times(boxprob_point(x-u,self.delta,self.h),boxprob_point(y-v,self.delta,self.h)),a))
        self.point_cache[key]=value
        return value
    def box_bound(self,i,B):
        xl,xr,yl,yr=B;h=self.h;d=self.delta
        px,pxx=cauchyjet(xl,xr,h);py,pyy=cauchyjet(yl,yr,h)
        xx=scale(times(pxx,py),self.lam);yy=scale(times(px,pyy),self.lam)
        for u,v,a in self.physical:
            tx,txx=boxjet(xl,xr,u,d,h);ty,tyy=boxjet(yl,yr,v,d,h)
            xx=plus(xx,scale(times(txx,ty),a));yy=plus(yy,scale(times(tx,tyy),a))
        # Sum SIGNED derivative intervals with ACTUAL positive coefficients,
        # then take max(0,sup). This preserves cancellations of curvature.
        Mx=max(0,xx[1]);My=max(0,yy[1])
        err=scale((Mx,Mx),(xr-xl)**2/F(8))[1]+scale((My,My),(yr-yl)**2/F(8))[1]
        corners=[(xl,yl),(xl,yr),(xr,yl),(xr,yr)]
        values=[minus(self.gpoint(x,y),iv(self.a[i])) for x,y in corners]
        lower=min(v[0] for v in values)-err
        pmin=times(cauchyjet(xl,xr,h)[0],cauchyjet(yl,yr,h)[0])[0]
        assert pmin>0
        repair_upper=max(F(0),F(-lower,pmin))
        neg=[(x,y,v[1]) for (x,y),v in zip(corners,values) if v[1]<0]
        return lower,err,repair_upper,neg,(xx,yy)

def split(B):
    xl,xr,yl,yr=B
    if xr-xl>=yr-yl:
        m=(xl+xr)/2
        return [(xl,m,yl,yr),(m,xr,yl,yr)]
    m=(yl+yr)/2
    return [(xl,xr,yl,m),(xl,xr,m,yr)]

def exact_run(args):
    start=time.monotonic();inp=Path(args.input)
    raw=inp.read_bytes();data=json.loads(raw,parse_float=str);repair=F(args.repair)
    assert repair>=0
    model=Model(data,repair);stack=[]
    for i,((x,y),_) in reversed(list(enumerate(model.oo))):
        for ix in range(args.initial_subdiv):
            for iy in range(args.initial_subdiv):
                xl=x-model.delta/2+model.delta*ix/args.initial_subdiv
                xr=xl+model.delta/args.initial_subdiv
                yl=y-model.delta/2+model.delta*iy/args.initial_subdiv
                yr=yl+model.delta/args.initial_subdiv
                stack.append((i,(xl,xr,yl,yr),0))
    tested=0;leaves=[];unresolved=[];witnesses=[];maxrepair=F(0);split_records=[]
    while stack and tested<args.max_boxes and time.monotonic()-start<args.max_seconds:
        i,B,depth=stack.pop();lo,err,req,neg,curv=model.box_bound(i,B);tested+=1
        record=dict(cell=i,box=list(map(str,B)),depth=depth,lower_fixed=str(lo),error_fixed=str(err),repair_upper=str(req),curvature_intervals_fixed=[[str(t)for t in a]for a in curv])
        if lo>=0:
            leaves.append(record)
            continue
        if neg:
            for x,y,upper in neg:
                # Closed-cell corner witness proves failure in its one-sided
                # interior neighborhood whenever strict upper<0 (continuity).
                witnesses.append(dict(cell=i,point=[str(x),str(y)],residual_upper_fixed=str(upper),one_sided_positive_area=True))
            if args.stop_on_negative:
                unresolved.append(record);maxrepair=max(maxrepair,req);break
        if depth>=args.max_depth:
            unresolved.append(record);maxrepair=max(maxrepair,req)
        else:
            split_records.append(dict(cell=i,box=list(map(str,B)),depth=depth))
            for child in split(B):stack.append((i,child,depth+1))
    complete=not stack
    status='PASS_SOURCE_NONNEGATIVE' if complete and not unresolved else ('INVALID_SOURCE' if witnesses else 'INCOMPLETE_SOURCE_CERTIFICATE')
    # Complete terminal coverage yields a sufficient EXTRA origin amplitude,
    # even when original source could not be certified. This is not D-certified.
    global_repair=str(maxrepair) if complete else None
    out=dict(status=status,input=str(inp.resolve()),input_sha256=hashlib.sha256(raw).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version,executable=sys.executable,bits=BITS,atan_terms=88,pi_fixed=list(map(str,PI)),rounding='Every signed arithmetic operation uses exact integers and outward floor/ceil at common scale 2^160; atan alternating remainder included.',h=str(model.h),delta=str(model.delta),N=model.N,origin_repair_added_exact=str(repair),source_origin_mass_exact=str(model.lam),mu_mass_exact=str(model.mass),halfmass_over_pi_exact=str(model.h*model.mass),source_only=True,D_minimum_certified=False,complete_terminal_coverage=complete,checked_boxes=tested,accepted_leaves=len(leaves),unresolved_leaves=len(unresolved),unchecked_stack=len(stack),additional_origin_repair_sufficient_exact=global_repair,elapsed_seconds=time.monotonic()-start,limits=dict(max_seconds=args.max_seconds,max_boxes=args.max_boxes,max_depth=args.max_depth,initial_subdiv=args.initial_subdiv),accepted=leaves,unresolved=unresolved,split_records=split_records,negative_corner_witnesses=witnesses,unchecked=[dict(cell=i,box=list(map(str,B)),depth=d)for i,B,d in stack],scope='Only positivity of lambda*p_h+P_h*f-f on all source cells; exterior positivity automatic. No D-minimum certificate or new upper constant.')
    dest=Path(args.output);assert str(dest.resolve()).startswith('/Volumes/data/')
    dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['status','input_sha256','origin_repair_added_exact','complete_terminal_coverage','checked_boxes','accepted_leaves','unresolved_leaves','unchecked_stack','additional_origin_repair_sufficient_exact','elapsed_seconds']},indent=2))
    return 0

def self_test(output):
    # Exact algebra and tiny whole-cell coverage, not an optimized candidate.
    assert times(iv(F(-2)),iv(F(3)))==iv(F(-6))
    assert divide(iv(F(-3)),iv(F(2)))==iv(F(-3,2))
    # f=0 gives g=lambda*p_h>0; 3x3 physical-cell support placeholders.
    d=dict(h='1',delta='1/8',N=1,source_origin_mass='1',density_coefficients=['0']*3)
    m=Model(d,F(0));checked=[]
    for i,((x,y),_) in enumerate(m.oo):
        B=(x-m.delta/2,x+m.delta/2,y-m.delta/2,y+m.delta/2)
        lo,err,req,neg,curv=m.box_bound(i,B)
        assert lo>0 and not neg and req==0
        checked.append(dict(cell=i,lower_fixed=str(lo),error_fixed=str(err)))
    # Nonzero density exact source via conservative point-density minorant.
    # 9 cells each area d²; choose a<=min p_h on support, so g-a>=p_h-a>=0.
    dd=F(1,8);aa=F(1,32);d['density_coefficients']=[str(aa)]*3
    m=Model(d,F(0));assert aa < 1/(F(22,7)**2*(1+F(3,16)**2)**2)
    nonzero=[]
    for i,((x,y),_) in enumerate(m.oo):
        B=(x-dd/2,x+dd/2,y-dd/2,y+dd/2)
        lo,err,req,neg,curv=m.box_bound(i,B)
        assert lo>0 and not neg
        nonzero.append(dict(cell=i,lower_fixed=str(lo),error_fixed=str(err)))
    out=dict(status='PASS_EXACT_TINY_SOURCE_SELF_TEST',bits=BITS,python=sys.version,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),pi_fixed=list(map(str,PI)),zero_density_boxes=checked,positive_density_boxes=nonzero,scope='Tiny exact rational source checks; no LP and no D lower bound.')
    p=Path(output);assert str(p.resolve()).startswith('/Volumes/data/');p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));return 0

def main():
    assert os.path.ismount('/Volumes/data')
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input');ap.add_argument('--output',required=True);ap.add_argument('--repair',default='0');ap.add_argument('--max-seconds',type=float,default=45);ap.add_argument('--max-boxes',type=int,default=20000);ap.add_argument('--max-depth',type=int,default=12);ap.add_argument('--initial-subdiv',type=int,default=1);ap.add_argument('--stop-on-negative',action='store_true');ap.add_argument('--self-test',action='store_true');args=ap.parse_args()
    assert args.initial_subdiv>=1 and args.max_depth>=0 and args.max_seconds>0 and args.max_boxes>0
    if args.self_test:return self_test(args.output)
    assert args.input
    return exact_run(args)
if __name__=='__main__':raise SystemExit(main())
