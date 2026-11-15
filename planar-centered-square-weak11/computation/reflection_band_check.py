#!/usr/bin/env python3
"""Exact Q(sqrt(2)) checks of the reflection-band witness, not a general proof."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
import json, platform, hashlib
from pathlib import Path

@dataclass(frozen=True)
class Q:
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,"a",F(self.a)); object.__setattr__(self,"b",F(self.b))
    def __add__(self,o):
        o=asq(o); return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,o): return self+-asq(o)
    def __rsub__(self,o): return asq(o)+-self
    def __mul__(self,o):
        o=asq(o); return Q(self.a*o.a+2*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=asq(o); den=o.a*o.a-2*o.b*o.b
        return self*Q(o.a/den,-o.b/den)
    def sign(self):
        a,b=self.a,self.b
        if not a:return (b>0)-(b<0)
        if not b:return (a>0)-(a<0)
        if a>0 and b>0:return 1
        if a<0 and b<0:return -1
        return ((a*a>2*b*b)-(a*a<2*b*b))*(1 if a>0 else -1)
    def __lt__(self,o):return (self-o).sign()<0
    def __le__(self,o):return (self-o).sign()<=0
    def __str__(self):return f"({self.a})+({self.b})sqrt(2)"
def asq(x): return x if isinstance(x,Q) else Q(x)
z=Q(); one=Q(1); sqrt2=Q(0,1); k=sqrt2-1; w=k/2
cases={
 "singleton":[(z,z)],
 "positive_gap_selected_away_from_c":[(z,z),(Q(F(1,10)),Q(-3))],
 "persistent_right_tie":[(z,z),(z,Q(3))],
 "distinct_top_bottom_radii":[(z,z),(Q(F(3,10)),Q(F(1,2)))],
 "just_inactive_pair":[(z,z),(z,2*sqrt2+Q(F(1,100)))],
 "activation_wall_empty_pair":[(z,z),(z,2*sqrt2)],
 "just_active_pair":[(z,z),(z,2*sqrt2-Q(F(1,100)))],
}
rows=[]
for name,pts in cases.items():
    rects=[]
    for size in range(1,len(pts)+1):
        rad=one if size==1 else sqrt2
        for idx in combinations(range(len(pts)),size):
            xs=[pts[j][0] for j in idx];ys=[pts[j][1] for j in idx]
            a,b,d,e=min(xs),max(xs),min(ys),max(ys)
            left,right,bottom,top=b-rad,a+rad,e-rad,d+rad
            if left<right and bottom<top:
                rects.append(dict(idx=idx,r=rad,a=a,b=b,d=d,e=e,left=left,right=right,bottom=bottom,top=top))
    c=max(R['right'] for R in rects)
    eligible=[R for R in rects if c-w<=R['right']]
    S=max(eligible,key=lambda R:R['top']);r=S['r'];Y=S['top'];g=c-S['right']
    Sminus=min(eligible,key=lambda R:R['bottom']);rminus=Sminus['r'];Z=Sminus['bottom'];gminus=c-Sminus['right']
    assert Z<Y
    assert z<=gminus<=w and one<=rminus and w<=k*rminus-gminus
    assert all(x<c for x,y in pts)
    assert z<=g<=w and one<=r
    assert w<=k*r-g
    # Mixed bbox extrema and its actual radius are calculated independently.
    mixed_left=min(S['a'],2*c-S['b']);mixed_right=max(S['b'],2*c-S['a'])
    radius=sqrt2*r
    assert mixed_right-radius==c+g-k*r
    assert mixed_left+radius==c-g+k*r
    assert S['d']+radius==Y+k*r
    assert S['e']-radius<Y-k*r
    for R in rects:
        assert R['right']<=c-w or (R['top']<=Y and Z<=R['bottom'])
    mixed_min=min(Sminus['a'],2*c-Sminus['b'])
    mixed_max=max(Sminus['b'],2*c-Sminus['a'])
    mixed_rad=sqrt2*rminus
    assert mixed_max-mixed_rad==c+gminus-k*rminus
    assert mixed_min+mixed_rad==c-gminus+k*rminus
    assert Sminus['e']-mixed_rad==Z-k*rminus
    assert Z+k*rminus<Sminus['d']+mixed_rad
    if name=='distinct_top_bottom_radii':
        assert S['idx']!=Sminus['idx'] and r!=rminus
    Rrad=r+rminus
    M=4*len(pts);gain=w*k*Rrad/M;floor=k*k*4/(4*M)
    assert floor<=gain
    D=16*sum((R['r'] for R in rects),z)+8*w+2*k*Rrad
    eta=w*k*Rrad/D
    assert z<eta<w/2 and eta<one/4
    assert 2*eta<k*r and 2*eta<k*rminus
    assert D*eta==w*k*Rrad
    assert eta==w*k*Rrad/(2*(8*sum((R['r'] for R in rects),z)+4*w+k*Rrad))
    retained=(2*w*k*Rrad-D*eta)/(2*M)
    assert retained==gain/2
    rows.append(dict(case=name,nonempty=len(rects),selected_top=list(S['idx']),selected_bottom=list(Sminus['idx']),c=str(c),gplus=str(g),gminus=str(gminus),rplus=str(r),rminus=str(rminus),gain=str(gain),eta=str(eta)))
assert F(71,25)**2>8
assert F(17,6)**2>8
assert (k*k-Q(F(4,25))).sign()>0
assert (k*k-Q(F(1,6))).sign()>0
d=F(1,12); rational_rows=[]
for N in range(1,5):
    assert d==F(2160,25920**(2**(N-1)))
    t=d/45;tau=d*d/180
    assert t<=F(1,5) and 9*t/(1-t)<=d/4 and tau/t==d/4
    next_d=d*d/2160
    assert next_d<=d/2
    rational_rows.append(dict(N=N,d=str(d),t=str(t),tau=str(tau)))
    d=next_d
geometric_rows=[]
for steps in range(1,9):
    geometric_sum=sum((F(1,2**j) for j in range(steps)),F(0))
    assert geometric_sum==2*(1-F(1,2**steps))
    geometric_rows.append(dict(steps=steps,sum=str(geometric_sum)))
out=dict(iterated_gain_geometric_sum=geometric_rows,rational_square_certificates={'71_over_25_squared_minus_8':str(F(71,25)**2-8),'17_over_6_squared_minus_8':str(F(17,6)**2-8)},rational_recursion=rational_rows,python=platform.python_version(),arithmetic='exact Fraction Q(sqrt(2))',threshold=1,atom_mass=4,random_seed=None,tolerance=None,cases=rows,status='PASS',scope='seven finite double-band witness checks; complete general proof is separate',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
print(json.dumps(out,indent=2))
