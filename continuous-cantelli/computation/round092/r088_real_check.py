"""Numerical calibration of constants and analytic moving-wall heat layer."""
import json
import math
import platform
import mpmath as mp
mp.mp.dps=55
p=lambda t,x:mp.exp(-x*x/(2*t))/mp.sqrt(2*mp.pi*t)
A,T,c=map(mp.mpf,('1.1','0.8','1.4'))
K=2*A/T*p(T,A)
barrier=[]
for lam in (25,100,400):
    est=lam*mp.quad(lambda v:mp.exp(-v)*(p(T,v/lam-A)-p(T,v/lam+A)),[0,2,8,mp.inf])
    barrier.append([lam,mp.nstr(est,19),mp.nstr(abs(est-K),8)])
edge=[]
for u in [mp.mpf(1)/20,mp.mpf(1)/100,mp.mpf(1)/500]:
    def difference(w):
        return (p(c*w,u-w)+p(c*w,-u-w)-2*p(c*w,-w))*p(1,w-A)
    val=-mp.quad(difference,[0,u*u,10*u*u,mp.mpf('.3')])/u
    edge.append([str(u),mp.nstr(val,20)])
C=mp.mpf('0.7')
far=[]
for lam in (40,160,640):
    def integrand(z):
        r=lam-z
        correction=A*lam*z/r +2*mp.log(lam) -(4*mp.log(r)+C)*lam**2/(2*r*r)
        return p(1,z)*mp.exp(correction)
    val=mp.quad(integrand,[-16,-4,0,4,16])
    far.append([lam,mp.nstr(val,22)])
drift=mp.mpf('-0.17')
walls=[]
for t in map(mp.mpf,('.73','.8','.87')):
    y=mp.mpf('-1.7')
    def layer(s):
        return (A/s)*p(s,A-drift*s)*p(t-s,y-(-A+drift*s))
    val=mp.quad(layer,[0,t/4,t/2,3*t/4,t])
    walls.append([str(t),mp.nstr(abs(val-p(t,y)),8)])
assert max(float(row[1]) for row in walls)<1e-45
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,dps=55,
                      random='none',A=str(A),T=str(T),c=str(c),C=str(C),
                      barrier_limit=mp.nstr(K,22),barrier=barrier,
                      near_edge_limit=mp.nstr(2*p(1,A)/c,22),near_edge=edge,
                      remote_limit=mp.nstr(mp.exp((A*A-C)/2),22),remote_window=[-16,16],remote=far,
                      moving_wall_exterior_errors=walls),indent=2))
