'''Rigorous positive-part dual integration on a parameter interval.'''
from fractions import Fraction as F
from flint import arb,arb_series,ctx
from certify_b1_classification import rat,dot,perp,unit

def normup(v):return sum((abs(x).upper()**2 for x in v),arb(0)).sqrt().upper()
def jets(G,w):
 t=arb_series([w,1,0]);a,b=G.ns['points'](t)[0][2];q=G.qa((a,b));p=(q[0][0],q[1][0]);v=(q[0][1],q[1][1]);acc=(2*q[0][2],2*q[1][2]);speed=dot(v,v).sqrt();assert speed>0
 u=tuple(-x/speed for x in perp(v));d=dot(v,acc);up=tuple(-x for x in perp(tuple(acc[i]/speed-v[i]*d/speed**3 for i in range(2))))
 return p,v,u,up

def calculate(G,raw,low,high):
 ctx.cap=3;pi=G.pi;mid=(low+high)/2;rad=(high-low)/2
 q,v,u,du=jets(G,rat(mid));qb,vb,ub,dub=jets(G,rat(low).union(rat(high)))
 dq=normup(vb)*rat(rad);dn=normup(dub)*rat(rad)
 q=tuple(x+arb(0,dq) for x in q);u=tuple(x+arb(0,dn) for x in u);direction=G.qa(u);height=dot(q,u)
 w0,wt,wf=map(rat,raw['w']);nu=list(map(rat,raw['nu']));assert w0>0 and wt>0 and wf>0
 th=61*pi/1800;den=(3*th).sin();sources=[];constant=wf/(2*den)-w0*G.height-wt*height;pc=arb(0)
 def add(vec,w):
  nonlocal constant,pc
  a=arb.atan2(vec[1],vec[0]);sgn=1
  if a<0:a+=2*pi
  if a>pi:a-=pi;sgn=-1
  assert 0<a<pi
  sources.append((a,w*sgn));constant+=w*(rat(F(1-sgn,2))+vec[0]/2);pc+=w*vec[1]
 add(G.normal,w0);add(direction,wt);add(unit(th),-wf/den);add(unit(th+pi/3),wf/den);add(unit(th+2*pi/3),-wf/den);add(unit(pi/3),-nu[0]);add(unit(2*pi/3),-nu[1]);constant+=nu[0]/2+nu[1]/2+nu[2]
 # Equality multiplier adjustment cancels the translation p exactly.
 shift=pc/(pi/3).sin();constant+=shift/4;sources.append((pi/3,-shift))
 edges=[arb(0)]+sorted((a for a,w in sources[:-1]),key=lambda x:float(x.mid()))+[pi]
 assert all(b>a for a,b in zip(edges,edges[1:]));total=arb(0)
 for l,r in zip(edges,edges[1:]):
  assert r-l<pi
  middle=(l+r)/2;aa=-nu[2];bb=-nu[3]
  for a,w in sources:
   assert a<middle or a>middle
   if a>middle:
    s,c=a.sin_cos();aa-=w*c;bb+=w*s
  kl=aa*l.sin()+bb*l.cos();kr=aa*r.sin()+bb*r.cos();Hl=-aa*l.cos()+bb*l.sin();Hr=-aa*r.cos()+bb*r.sin();radius=normup((aa,bb));options=[]
  if not kl<0 and not kr<0:options.append(Hr-Hl)
  if not kl>0 and not kr>0:options.append(arb(0))
  if not kl<0 and not kr>0:options.append(radius-Hl)
  if not kl>0 and not kr<0:options.append(Hr+radius)
  assert options;total+=max(x.upper() for x in options)
 return constant+total
