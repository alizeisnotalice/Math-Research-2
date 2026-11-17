"""Rational interval certification of a fixed-angle candidate path in R."""
from fractions import Fraction as F
from math import isqrt
import json,sys,hashlib
from pathlib import Path
from rational_trig import PI,trig,scale
Q=10**24
class I:
 def __init__(self,a,b=None):
  if isinstance(a,I):self.lo,self.hi=a.lo,a.hi;return
  self.lo=F(a);self.hi=F(a if b is None else b)
  assert self.lo<=self.hi
 def rounded(self):return I(F((self.lo*Q).__floor__(),Q),F((self.hi*Q).__ceil__(),Q))
 def __add__(self,b):
  b=I(b);return I(self.lo+b.lo,self.hi+b.hi).rounded()
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,b):return self+-I(b)
 def __rsub__(self,b):return I(b)+-self
 def __mul__(self,b):
  b=I(b);v=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)];return I(min(v),max(v)).rounded()
 __rmul__=__mul__
 def __truediv__(self,b):
  b=I(b);assert b.lo*b.hi>0
  return self*I(1/b.hi,1/b.lo).rounded()
 def sqrt(self):
  assert self.lo>=0
  lo=isqrt((self.lo*Q*Q).__floor__());hi=isqrt((self.hi*Q*Q).__floor__())+1
  assert F(lo*lo,Q*Q)<=self.lo and F(hi*hi,Q*Q)>=self.hi
  return I(F(lo,Q),F(hi,Q))
 def out(self):return [str(self.lo),str(self.hi)]
def V(x,y):return (I(x),I(y))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def mul(t,a):return tuple(t*x for x in a)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),I(0))
def perp(a):return (-a[1],a[0])
def norm2(a):return dot(a,a)
def upper_cross(l,n):
 q=sub(n,l);d=norm2(q).sqrt();h=(1-norm2(q)/4).sqrt()
 return sub(mul(I(F(1,2)),add(l,n)),mul(h/d,perp(q)))
r3=I(3).sqrt();half=I(F(1,2));a=scale(PI,F(5,6)+F('1.5494')/180)
co=I(*trig(a,'cos'));si=I(*trig(a,'sin'))
F3=V(-half,(half+half*co)/si)
E=V(half*F3[0]-r3*half*F3[1],r3*half*F3[0]+half*F3[1])
C=V(-E[0],E[1]);n=V(half,-r3*half)
projection=dot(n,F3)-half
foot=sub(F3,mul(projection,n))
G=add(foot,mul((1-projection*projection).sqrt(),V(r3*half,half)))
L0=V(r3*half-half,((r3*half-half)-1)/r3)
M=V(F(-1,4),-r3/4);Nm=sub(mul(I(2),M),E)
Z=upper_cross(L0,M)
Y=V(0,C[1]+(1-C[0]*C[0]).sqrt())
z=sub(Z,C);d2=norm2(z)
assert d2.lo>1
T=add(C,sub(mul(I(1)/d2,z),mul((d2-1).sqrt()/d2,perp(z))))
assert T[0].hi<0 and Z[0].hi<T[0].lo
assert T[1].lo>C[1].hi
# Basic angle range and ordering for L(s), established without inverse trig.
v=sub(G,F3);e=V(r3*half,-half)
assert v[1].hi<0 and v[0].lo>0
assert (e[0]*v[1]-e[1]*v[0]).lo>0
assert C[0].hi<G[0].lo
# Endpoint disk membership. T is exactly on C, and tangent segment stays outside C.
assert norm2(sub(Z,G)).hi<1 and norm2(sub(T,G)).hi<1
assert norm2(sub(Z,E)).hi<1
# Arc uses x from Tx to zero, y=Cy+sqrt(1-(x-Cx)^2).
arc_bounds=[]
for j in range(64):
 u=I(F(j,64),F(j+1,64));x=T[0]*u
 y=C[1]+(1-(x-C[0])*(x-C[0])).sqrt()
 b=norm2(sub(V(x,y),G));assert b.hi<1
 arc_bounds.append(b.hi)
# Arc and segment are inside this coarse rectangle; endpoint and arc monotonicity suffice.
box=V(I(F(-45,1000),0),I(F(545,1000),F(564,1000)))
for P in (Y,Z,T):
 for x,b in zip(P,box):assert b.lo<x.lo or b.lo==x.lo;assert x.hi<=b.hi
# Circle arc y increases with x on this interval, so endpoints bound all y.
assert T[0].hi<C[0].lo
assert norm2(sub(box,Nm)).lo>1
# Unique inverse centers on segments C--G and Nm--E: distances strictly decrease along both.
Lbox=tuple(I(min(c.lo,g.lo),max(c.hi,g.hi)) for c,g in zip(C,G))
Nbox=tuple(I(min(c.lo,g.lo),max(c.hi,g.hi)) for c,g in zip(Nm,E))
assert norm2(sub(Lbox,Nbox)).lo>0 and norm2(sub(Lbox,Nbox)).hi<4
dL=V(r3*half,half);dN=V(-r3*half,half)
assert dot(sub(box,Lbox),dL).lo>0
assert dot(sub(box,Nbox),dN).lo>0
assert box[1].lo>max(Lbox[1].hi,Nbox[1].hi)
out={'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'trig_helper_sha256':hashlib.sha256(Path(__file__).with_name('rational_trig.py').read_bytes()).hexdigest(),'sigma_degrees_exact':'1.5494','rounding_denominator':Q,'arc_cells':64,'coordinates':{name:[x.out() for x in p] for name,p in [('F3',F3),('E3',E),('C',C),('G',G),('Nm',Nm),('Y',Y),('Z',Z),('T',T)]},'arc_G_distance_squared_upper':str(max(arc_bounds)),'whole_path_Nm_distance_squared_lower':str(norm2(sub(box,Nm)).lo),'inverse_L_monotonicity_lower':str(dot(sub(box,Lbox),dL).lo),'inverse_N_monotonicity_lower':str(dot(sub(box,Nbox),dN).lo),'passed':True,'scope':'continuous candidate path in parameter image R at one fixed sigma; not a universal-cover or area certificate'}
print(json.dumps(out,indent=2))
