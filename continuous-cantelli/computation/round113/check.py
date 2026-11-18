"""Deterministic calibration of the cutoff and an oscillatory sharpness caution."""
import json
import math
import platform
import numpy as np
import scipy
import sympy as sp
from scipy.integrate import quad
from scipy.special import i0

u=sp.symbols('u')
fourier=sp.expand(((8*u+4)**3-4*(8*u+2)**3+6*(8*u)**3)/32)
assert -sp.diff(fourier,u,2).subs(u,0)==48

def khat(s):
    z=abs(8*s)
    if z<=2:
        return 1-3*z*z/8+3*z**3/32
    if z<4:
        return (4-z)**3/32
    return 0.

def density(t,y):
    return quad(lambda v: math.cos(v*y)*math.exp(-v*v/2)*khat(v/t)/math.pi,
                0,min(16.,t/2),epsabs=2e-13,epsrel=2e-13)[0]

b,k=.4,.7
a=2*math.log(i0(b/2))
def profile(x):
    return 1+(a+b*math.sin(k*x*x))/(1+x*x)
rows=[]
for t in (10,50,100,300):
    length=quad(lambda x:math.exp(-x*x*(profile(x)-1)/2),t,t+1,
                epsabs=1e-11,epsrel=1e-11,limit=1000)[0]
    rows.append(dict(t=t,window_average=length))
assert abs(rows[-1]['window_average']-1)<.002
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,numpy=np.__version__,
                      sympy=sp.__version__,deterministic=True,cutoff='normalized sinc^4(x/8)',
                      near_origin_fourier=str(fourier),noise_second_moment=48,
                      gaussian_cutoff_density=[dict(t=t,y0=density(t,0),y1=density(t,1)) for t in (4,16,64)],
                      auxiliary_profile='not asserted exact: 1+(a+b sin(k x^2))/(1+x^2)',
                      a=a,b=b,k=k,liminf_f_minus_x_fprime=1-2*k*b,rows=rows),indent=2))
