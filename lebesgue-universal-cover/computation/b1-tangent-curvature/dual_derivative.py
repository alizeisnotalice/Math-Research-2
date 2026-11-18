'''Interval first derivative of the explicitly integrated continuous dual.'''
from fractions import Fraction as F
from flint import arb,arb_series,ctx
from certify_b1_classification import rat,dot,perp,unit

def constant(x):return arb_series([x])
def derivative_bound(G,raw,low,high):
 ctx.cap=3;pi=G.pi;W=arb_series([rat(low).union(rat(high)),1]);a,b=G.ns['points'](W)[0][2];q=G.qa((a,b));v=tuple(x.derivative() for x in q);speed=dot(v,v).sqrt();u=tuple(-x/speed for x in perp(v));direction=G.qa(u);height=dot(q,u)
 w0,wt,wf=map(rat,raw['w']);nu=list(map(rat,raw['nu']));th=61*pi/1800;den=(3*th).sin();sources=[];C=constant(wf/(2*den)-w0*G.height)-wt*height;pc=constant(arb(0))
 def add(vec,weight):
  nonlocal C,pc
  x,y=vec;x=constant(x) if not isinstance(x,arb_series) else x;y=constant(y) if not isinstance(y,arb_series) else y
  a=arb.atan2(y[0],x[0]);sgn=1
  if a<0:a+=2*pi
  if a>pi:a-=pi;sgn=-1
  assert 0<a<pi
  ad=(x[0]*y[1]-y[0]*x[1])/(x[0]**2+y[0]**2)
  angle=arb_series([a,ad]);sources.append((angle,weight*sgn));C+=weight*(rat(F(1-sgn,2))+x/2);pc+=weight*y
 add(G.normal,w0);add(direction,wt);add(unit(th),-wf/den);add(unit(th+pi/3),wf/den);add(unit(th+2*pi/3),-wf/den);add(unit(pi/3),-nu[0]);add(unit(2*pi/3),-nu[1]);C+=nu[0]/2+nu[1]/2+nu[2]
 shift=pc/(pi/3).sin();C+=shift/4;sources.append((constant(pi/3),-shift))
 edges=[constant(arb(0))]+sorted((a for a,w in sources[:-1]),key=lambda a:float(a[0].mid()))+[constant(pi)]
 assert all(b[0]>a[0] for a,b in zip(edges,edges[1:])),[(str(a[0]),str(b[0])) for a,b in zip(edges,edges[1:]) if not b[0]>a[0]]
 deriv=C[1]
 for l,r in zip(edges,edges[1:]):
  assert r[0]-l[0]<pi
  mid=(l[0]+r[0])/2;aa=constant(-nu[2]);bb=constant(-nu[3])
  for a,w in sources:
   assert a[0]<mid or a[0]>mid
   if a[0]>mid:sn,cs=a.sin_cos();aa-=w*cs;bb+=w*sn
  sl,cl=l.sin_cos();sr,cr=r.sin_cos();kl=aa*sl+bb*cl;kr=aa*sr+bb*cr;Hl=-aa*cl+bb*sl;Hr=-aa*cr+bb*sr;opt=[]
  if not kl[0]<0 and not kr[0]<0:opt.append((Hr-Hl)[1])
  if not kl[0]>0 and not kr[0]>0:opt.append(arb(0))
  if (not kl[0]<0 and not kr[0]>0) or (not kl[0]>0 and not kr[0]<0):
   radius=(aa*aa+bb*bb).sqrt()
   if not kl[0]<0 and not kr[0]>0:opt.append((radius-Hl)[1])
   if not kl[0]>0 and not kr[0]<0:opt.append((Hr+radius)[1])
  assert opt and all(x.is_finite() for x in opt)
  v=opt[0]
  for x in opt[1:]:v=v.union(x)
  deriv+=v
 return deriv
