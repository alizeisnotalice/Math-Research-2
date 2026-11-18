"""Rational regression checks for the exact grouped area identity and perturbation enclosure."""
from random import Random
from fractions import Fraction as F
from collections import defaultdict
rng=Random(880103)
# Exact rational complex arithmetic (no binary complex coercions).
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def rot(a,z):return (a[0]*z[0]-a[1]*z[1],a[0]*z[1]+a[1]*z[0])
for trial in range(500):
 n=rng.randint(3,20);owners=rng.choices(range(5),k=n);p=[(rng.randint(-7,7),rng.randint(-7,7)) for i in range(n)]
 t=[(rng.randint(-7,7),rng.randint(-7,7)) for k in range(5)];u=[(rng.randint(-6,6),rng.randint(-6,6)) for k in range(5)]
 z=[];theta=[]
 for k in range(5):
  x=F(rng.randint(-10,10),11);z.append(((1-x*x)/(1+x*x),2*x/(1+x*x)));theta.append(2*abs(x))
 base=[add(t[k],a) for k,a in zip(owners,p)];act=[add(add(t[k],u[k]),rot(a,z[k])) for k,a in zip(owners,p)]
 def twice(q):return sum(det(a,b) for a,b in zip(q,q[1:]+q[:1]))
 exact=twice(act)-twice(base)
 g=defaultdict(lambda:(0,0));gr=defaultdict(int);co=defaultdict(int);mix=defaultdict(lambda:(0,0));cnt=defaultdict(int);cr=defaultdict(int);dt=defaultdict(int)
 for i in range(n):
  c=(i+1)%n;k=owners[i];j=owners[c];a=p[i];b=p[c];d=sub(base[c],base[i-1])
  g[k]=add(g[k],d);gr[k]-=dot(a,d);co[k]+=det(a,d)
  mix[k,j]=add(mix[k,j],b);mix[j,k]=sub(mix[j,k],a)
  if k==j:co[k]-=2*det(a,b)
  else:
   f,l=sorted([k,j]);sgn=1 if k<j else -1
   cnt[f,l]+=sgn;cr[f,l]+=det(a,b);dt[f,l]+=sgn*dot(a,b)
 calc=0;error=0
 for k in range(5):
  ck=z[k][0]-1;sk=z[k][1]
  calc+=det(u[k],g[k])+gr[k]*sk+co[k]*ck
  hx,hy=map(abs,u[k]);error+=hx*abs(g[k][1])+hy*abs(g[k][0])+abs(gr[k])*theta[k]+max(co[k],0)*theta[k]**2/2
  for j in range(5):
   v=mix[k,j];calc+=det(u[k],sub(rot(v,z[j]),v))
   error+=theta[j]*(hx*abs(v[0])+hy*abs(v[1]))+theta[j]**2/2*(hx*abs(v[1])+hy*abs(v[0]))
  for j in range(k+1,5):
   cj=z[j][0]-1;sj=z[j][1]
   calc+=cnt[k,j]*det(u[k],u[j])+cr[k,j]*(ck*cj+sk*sj)+dt[k,j]*(ck*sj-sk*cj)
   error+=abs(cnt[k,j])*(hx*abs(u[j][1])+hy*abs(u[j][0]))
   error+=theta[k]*theta[j]*(abs(cr[k,j])+abs(dt[k,j])*min(1,(theta[k]+theta[j])/2))
 assert calc==exact and exact>=-error,(trial,calc,exact,error)
print('500 exact polygon identities and 500 one-sided perturbation enclosures passed (seed 880103)')
