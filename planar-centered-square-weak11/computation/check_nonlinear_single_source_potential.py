#!/usr/bin/env python3
"""Exact checks for nonlinear single-source tails and the hard-core pair certificate.
Q2 arithmetic/slab union reuse the independently replayed R32 implementation;
new geometry and potential checks below do not modify the frozen R32 files.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
from itertools import product, combinations
from pathlib import Path
import hashlib, json, os, platform, sys

@total_ordering
@dataclass(frozen=True)
class Q2:
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))
    @staticmethod
    def coerce(x):
        return x if isinstance(x,Q2) else Q2(F(x))
    def __add__(self,x):
        x=self.coerce(x); return Q2(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return Q2(-self.a,-self.b)
    def __sub__(self,x): return self+-self.coerce(x)
    def __rsub__(self,x): return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x); return Q2(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self.coerce(x)
        den=x.a*x.a-2*x.b*x.b
        if den==0: raise ZeroDivisionError
        return self*Q2(x.a/den,-x.b/den)
    def sign(self):
        a,b=self.a,self.b
        if not b: return (a>0)-(a<0)
        if not a: return (b>0)-(b<0)
        if a>0 and b>0: return 1
        if a<0 and b<0: return -1
        s=(a*a>2*b*b)-(a*a<2*b*b)
        return s if a>0 else -s
    def __lt__(self,x): return (self-self.coerce(x)).sign()<0
    def __eq__(self,x):
        x=self.coerce(x); return self.a==x.a and self.b==x.b
    def record(self):
        return {'rational':str(self.a),'sqrt2_coefficient':str(self.b)}

ZERO=Q2(); ONE=Q2(1); HALF=Q2(F(1,2)); G=Q2(3,-2); OPT=Q2(4,-2)
assert G>0

def rectangle_union(rects):
    """Exact slab union; clipped endpoints stay in Q(sqrt(2))."""
    xs=sorted({x for r in rects for x in r[:2]})
    area=ZERO
    for left,right in zip(xs,xs[1:]):
        middle=(left+right)/2
        intervals=sorted((lo,hi) for x0,x1,lo,hi in rects if x0<middle<x1)
        total=ZERO
        if intervals:
            lo,hi=intervals[0]
            for a,b in intervals[1:]:
                if a<=hi: hi=max(hi,b)
                else: total+=hi-lo; lo,hi=a,b
            total+=hi-lo
        area+=(right-left)*total
    return area


M=F(9,4)
B=Q2(F(9,4),F(-3,2))
Z=Q2(3,F(-4,3))
assert B>Q2(F(1,8))

def a_rational(j,t,h):
    return 2*t*max(F(0),min(h,j*h+t)-max(F(0),j*h-t))

def phi_integer_radius(j,l,h):
    return 2*l*h*h*max(0,l-j) if j>=0 else -2*l*h*h*max(0,l+j)

def full_hardcore_strip(window,h=F(1),strip=True):
    # True one-row finite geometry. T=8h and R=2h isolate every vertical lift.
    pts=[j*h for j,c in zip(range(-3,4),window) if c]
    rects=[];boxes=0;surviving=0
    for left_index,left in enumerate(pts):
      for right_index in range(left_index,len(pts)):
        right=pts[right_index];count=right_index-left_index+1
        mass=M*h*h*count;span=right-left;boxes+=1
        # Exact width/non-cap feasibility before taking any root.
        if span*span>=mass or span>=4*h:continue
        assert count in (1,2)
        if count==2:assert span==2*h
        radius=Q2(F(3,4)*h) if count==1 else Q2(0,F(3,4)*h)
        a=Q2(right)-radius;b=Q2(left)+radius
        if strip:a=max(Q2(0),a);b=min(Q2(h),b)
        if a<b:
            rects.append((a,b,-radius,radius));surviving+=1
    return rectangle_union(rects),boxes,surviving

def main():
    assert os.path.ismount('/Volumes/data')
    harmonic_checks=0;h=F(1,64)
    for l in range(1,33):
      m=4*l*l*h*h;t=l*h
      for j in range(-65,66):
        actual=a_rational(j,t,h)
        assert actual==2*l*h*h*(1-l<=j<=l)
        diff=phi_integer_radius(j-1,l,h)-phi_integer_radius(j,l,h)
        assert diff==actual-m*(j==0)
        assert abs(phi_integer_radius(j,l,h))<=m/2
        harmonic_checks+=1
    # A cap-active source costs its true mass, with exactly unused m-A at j=0.
    h=F(1);t=F(2);A=4*t*t;m=F(25)
    def cap_phi(j):
        return sum((a_rational(q,t,h) for q in range(j+1,6)),F(0))-A*(j<0)
    for j in range(-3,4):
        assert cap_phi(j-1)-cap_phi(j)==a_rational(j,t,h)-A*(j==0)
        assert cap_phi(j-1)-cap_phi(j)-(a_rational(j,t,h)-m*(j==0))==(m-A)*(j==0)
    periodic_tail_levels=[];periodic_tail_checks=0
    for n in range(1,10):
        T=F(2**n);R=T;maximum_mass=T;quantum=F(1,2**(4*n))
        assert maximum_mass<=T*T/2
        assert maximum_mass/4<T*T/4 and maximum_mass/4<R*R
        assert (maximum_mass/quantum).denominator==1
        tail_at_T_over_2=2*maximum_mass/(T*T)
        assert tail_at_T_over_2<=1 and tail_at_T_over_2==F(1,2**(n-1))
        for ratio in (F(1,2),F(3,4),F(1),F(3,2),F(2),F(4)):
            r=ratio*T
            bound=maximum_mass/(2*T*r)+maximum_mass/(4*r*r)
            assert bound<=tail_at_T_over_2<=1
            periodic_tail_checks+=1
        periodic_tail_levels.append({'n':n,'T':str(T),'R':str(R),
          'maximum_mass':str(maximum_mass),'allowed_upper_bound_T_squared_over_2':str(T*T/2),
          'tail_at_T_over_2':str(tail_at_T_over_2),'cap_inactive_at_maximum_mass':True})
    legal=lambda v:all(not(a and b) for a,b in zip(v,v[1:]))
    vertices=[v for v in product((0,1),repeat=6) if legal(v)]
    windows=[v for v in product((0,1),repeat=7) if legal(v)]
    assert len(vertices)==21 and len(windows)==34
    records=[];geometry_evaluations=0;box_evaluations=0
    for scale in (F(1,8),F(1),F(7,3)):
      for w in windows:
        actual,boxes,surviving=full_hardcore_strip(w,scale)
        expected=(Q2(M/2)*(w[3]+w[4])+B*(w[2]*w[4]+w[3]*w[5]))*(scale*scale)
        assert actual==expected,(scale,w,actual,expected)
        pu=(Q2(M/2)*w[3]+B*w[1]*w[3])*(scale*scale)
        pv=(Q2(M/2)*w[4]+B*w[2]*w[4])*(scale*scale)
        slack=pv-pu-(actual-Z*(M*scale*scale*w[3]))
        assert slack==B*(scale*scale)*w[3]*(2-w[1]-w[5]) and slack>=0
        assert w[:-1] in vertices and w[1:] in vertices
        if scale==1:records.append({'window':w,'area':actual.record(),'potential_slack':slack.record(),
          'source_boxes':boxes,'receiving_rectangles':surviving})
        geometry_evaluations+=1;box_evaluations+=boxes
    pair=(0,0,1,0,1,0,0)  # positions -1 and +1, center empty
    local,_,_=full_hardcore_strip(pair)
    glob,_,_=full_hardcore_strip(pair,strip=False)
    assert local==Q2(M/2)+B
    assert glob==Q2(2*M)+2*B
    # Explicit open U=(0,1/32)x(-1/32,1/32), one radius for all centers.
    r=F(33,32)
    assert 2*M>4*r*r and F(31,32)>F(3,4)
    assert (2*M)/(4*r*r)==F(128,121)
    assert F(1,32)*F(1,16)==F(1,512)
    period_windows=[(0,1,0,1,0,1,0),(1,0,1,0,1,0,1)]
    period_area=sum((full_hardcore_strip(w)[0] for w in period_windows),ZERO)
    period_mass=M*sum(w[3] for w in period_windows)
    assert period_area/period_mass==Z
    assert F(9,64)+F(45,128)+F(9,64)==F(81,128)<1
    source=Path(__file__).resolve()
    result={'status':'PASS','python_version':sys.version,'platform':platform.platform(),
      'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'arithmetic':'Fraction and exact Q(sqrt(2)); no floating comparison','seed':None,'tolerance':None,
      'single_source_identity_checks':harmonic_checks,'cap_active_checks':7,
      'R32_periodic_tail_checks':periodic_tail_checks,'R32_periodic_tail_levels':periodic_tail_levels,
      'hardcore_vertices':len(vertices),'hardcore_edges':len(windows),
      'physical_scales':['1/8','1','7/3'],'geometry_evaluations':geometry_evaluations,
      'source_boxes_evaluated':box_evaluations,'B_at_h1':B.record(),'z_hardcore':Z.record(),
      'empty_center_pair_strip_area':local.record(),'pair_global_area':glob.record(),
      'strict_witness_area':'1/512','strict_witness_radius':'33/32','strict_witness_average':'128/121',
      'restricted_full_scale_tail':'81/128','base_scale_windows':records,
      'scope':'all single-source identities in the stated finite test family and complete hard-core local geometry; universal statements use analytic proof; not a general C2 upper bound'}
    out=source.parent/'results/nonlinear-single-source-potential/exact-check.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','single_source_identity_checks','cap_active_checks','R32_periodic_tail_checks',
      'hardcore_vertices','hardcore_edges','geometry_evaluations','source_boxes_evaluated',
      'B_at_h1','z_hardcore','empty_center_pair_strip_area','pair_global_area',
      'strict_witness_area','strict_witness_average','restricted_full_scale_tail']},indent=2))
    print(out)
if __name__=='__main__':main()
