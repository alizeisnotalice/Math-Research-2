"""Finite numerical check of the uniform localized multiplier error in (4)."""
import json
import math
import platform
import numpy as np
import scipy
from scipy.integrate import quad

def cutoff(u):
    a=abs(u)
    if a<=.25: return 1.0
    if a>=.5: return 0.0
    q=4*(a-.25)
    eta=lambda x: math.exp(-1/x) if x>0 else 0.
    return eta(1-q)/(eta(q)+eta(1-q))

def kernel(t,A,z):
    def integrand(s):
        scaled=math.exp(-s*s/2)/math.sqrt(2*math.pi)
        return scaled*cutoff(s/t)*math.exp(-A*(1+s/t)**2)
    edges=[-t*.5,-t*.25,0,t*.25,t*.5]
    # The reference bump φ(z)=exp(-z²/2) has Gaussian Fourier amplitude.
    re=sum(quad(lambda s:integrand(s)*math.cos(s*z),x,y,epsabs=1e-14,epsrel=1e-13,limit=200)[0] for x,y in zip(edges,edges[1:]))
    im=sum(quad(lambda s:integrand(s)*math.sin(s*z),x,y,epsabs=1e-14,epsrel=1e-13,limit=200)[0] for x,y in zip(edges,edges[1:]))
    err=abs(complex(re,im)-math.exp(-A-z*z/2))
    return t*(1+z*z)*math.exp(A/16)*err

tvals=[2,4,8,16,32,64]
avals=[0,.01,.25,1,4,16,64]
zvals=[0,.5,1.,2.,4.,6.]
report=[dict(t=t,grid_max_normalized_error=max(kernel(t,A,z) for A in avals for z in zvals)) for t in tvals]
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,A_values=avals,z_values=zvals,frequency_rows=report,formula='t*(1+z^2)*exp(A/16)*abs(error)',quadrature_abs_tol=1e-14,quadrature_rel_tol=1e-13),indent=2))
