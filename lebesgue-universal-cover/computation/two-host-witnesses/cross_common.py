"""Shared, pinned analytic supports and fixed hosts for the cross certificate."""
if not __debug__:
 raise RuntimeError("Run without Python -O: certificate checks use assertions")
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import sys,json,time,platform
import flint
from flint import arb,ctx
# Set precision before importing source modules that create Arb constants.
ctx.prec=160
C=Path(sys.argv[1]).resolve();D=Path(__file__).resolve().parent
# Pin all source imports and dynamically referenced source files to this replay.
manifest_path=D/'source-manifest.json'
assert sha256(manifest_path.read_bytes()).hexdigest() == '692121cc350edf4b3b2e3bf750d22e2a704fcd1210f2ce725405a82a14e27839'
manifest=json.loads(manifest_path.read_text())
for name, digest in manifest.items():
 assert sha256((C/name).read_bytes()).hexdigest()==digest, name
sys.path.insert(0,str(C))
from certify_multi_host_points import SupportHull,q
from certify_b9_padding_four_minimizer import axis,balance
from deformed_reuleaux_support import body,environment
assert ctx.prec == 160
pi=arb.pi(); A=lambda s:q(F(str(s)))
def dot(u,v):return u[0]*v[0]+u[1]*v[1]
def rot(a,v):return (a.cos()*v[0]-a.sin()*v[1],a.sin()*v[0]+a.cos()*v[1])
def neg(v):return (-v[0],-v[1])
ns=environment(str(C),'1.5497'); assert ctx.prec == 160
_,hp=body(ns,F('.7982907'),F('.032688')); assert ctx.prec == 160
# Opposite sectors use opposite edge directions and the same vertex with
# complementary arc flags: their support values sum exactly to one.
assert len(ns['orders']) == len(ns['sector']) == 10
assert all(ns['orders'][i+5] == ns['orders'][i][::-1] for i in range(5))
assert all(ns['sector'][i][0] == ns['sector'][i+5][0] and ns['sector'][i][1]+ns['sector'][i+5][1] == 1 for i in range(5))
L=pi/3;p=A('.1391310926383');x2=A('626.5943632307')*pi/2880;x1=L-((L-x2).cos()-arb(1)/2+L.cos()/2-p*L.sin()).acos();x5=2*L+(arb(3)/4+p*L.sin()).acos();unit=lambda t:(t.cos(),t.sin())
def sub(x,y):return(x[0]-y[0],x[1]-y[1])
def mul(x,k):return(x[0]*k,x[1]*k)
def add(x,y):return(x[0]+y[0],x[1]+y[1])
D0=sub(sub(unit(x2),unit(x1)),unit(x5));m=arb.atan2(D0[0],-D0[1]);ell=(dot(D0,D0).sqrt()/2).asin();xs=[arb(0),x1,x2,m-ell,m+ell,x5,pi];assert all(x<y for x,y in zip(xs,xs[1:]));assert x2<L<xs[3]<xs[4]<2*L<x5;point=(arb(1)/2,p);pieces=[]
for i,(a,b) in enumerate(zip(xs[:-1],xs[1:])):
 r=1-i%2;c=sub(point,mul(unit(a),r));pieces.append((r,c));point=add(c,mul(unit(b),r))
assert all(z.contains(0) for z in sub(point,sub((arb(1)/2,p),unit(arb(0)))))
def hf(v):
 t=arb.atan2(v[1],v[0])
 if t>=0: inv=False
 elif t<0: inv=True; t=t+pi; v=neg(v)
 else: raise ValueError(('angle',t))
 cand=[r+dot(c,v) for (a,b),(r,c) in zip(zip(xs[:-1],xs[1:]),pieces) if not (t<a or t>b)]
 assert cand
 h=cand[0]
 for z in cand[1:]: h=h.union(z)
 return 1-h if inv else h
# Six support bounds in one slightly rotated regular hexagon prove |x| < .58.
# The shift avoids evaluating directions exactly at a switching angle.
eps=A('0.000000000001')
for name,fn in [('pent',hp),('five',hf)]:
 for i in range(6):
  h=fn(unit(pi*i/3+eps))
  assert h<A('.50000000001'),(name,i,h)
assert 2*A('.50000000001')/arb(3).sqrt()<A('.58')
assert sha256((D/'certificate-labels.json').read_bytes()).hexdigest() == '5249239de976c2380eb8d3dfd3bb47347ab3cbf48ca11f3188fb07121a5a4da5'
class Host:
 def __init__(self,pose,seven):
  self.h=SupportHull(list(map(lambda x:q(F(x)),pose[:11])));self.points=[]
  for name,key,j in [('finite-point-leaf-input.json','denominator',11),('multi-host-points.json','point_denominator',14)]:
   if j==14 and not seven:continue
   d=json.loads((C/name).read_text());points=[(q(F(x,d[key])),q(F(y,d[key]))) for x,y in d['points']]
   if j==14:
    mean=tuple(sum(z[k] for z in points)/len(points) for k in range(2));points=[((z[0]-mean[0]),q(F(pose[17]))*(z[1]-mean[1])) for z in points]
   a=q(F(pose[j]));shift=(q(F(pose[j+1])),q(F(pose[j+2])));self.points.extend(add(rot(a,z),shift) for z in points)
 def row(self,seed):
  w,u,h=self.h.direction(seed)
  for z in self.points:h=max(h,(dot(u,z)).upper())
  return w,u,h
