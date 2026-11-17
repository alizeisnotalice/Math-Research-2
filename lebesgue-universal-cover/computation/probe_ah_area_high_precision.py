"""Independent high-precision adaptive integration, numerical cross-check only."""
import runpy,contextlib,io,json,sys
from pathlib import Path
import mpmath as m
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(Path(__file__).with_name('check_ah_high_precision.py')))
m.mp.dps=60
C,G,E,M,Y,Z=(d[k] for k in ('C','G','E3','M','Y','Z'))
h=m.mpf(1)/2;s3=m.sqrt(3);J=lambda v:m.matrix([-v[1],v[0]])
nL=m.matrix([h,-s3/2]);nN=m.matrix([-h,-s3/2]);eL=m.matrix([s3/2,h]);eN=m.matrix([-s3/2,h])
z=Z-C;d2=m.fdot(z,z);T=C+z/d2-m.sqrt(d2-1)*J(z)/d2
upper=d['X']
def pair(x):
 def inv(n,e):
  q=m.fdot(n,x)-h
  return x-q*n-m.sqrt(1-q*q)*e
 L=inv(nL,eL);N=inv(nN,eN);c=L[0]+h
 lx=c/2+s3/2*m.sqrt(1-c*c)-h
 return upper(m.matrix([lx,(lx-1)/s3]),2*M-N)
def line(w):return pair(Z+w*(T-Z))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def density(w):
 p=line(w);dp=m.matrix([m.diff(lambda t:line(t)[j],w) for j in (0,1)])
 return cross(p,dp)
def arc(a,b,c):
 u=a-c;v=b-c
 return cross(c,b-a)+m.atan2(cross(u,v),m.fdot(u,v))
Tp=pair(T);Xas=upper(G,E)
fixed=arc(Tp,Xas,G)+arc(Xas,Y,E)+arc(Y,T,C)+cross(T,Z)
area=abs((fixed+m.quad(density,[0,m.mpf('.5'),1]))/2)
print(json.dumps({'python':sys.version,'mpmath_version':m.__version__,'precision':60,'AH_area_numerical':m.nstr(area,55),'scope':'adaptive high precision cross-check only; rational interval certificate carries rigor'},indent=2))
