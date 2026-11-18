"""Finite Brownian-grid and half-line-kernel calibration; no exact graph is assumed."""
import json
import math
import platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import ndtr
from numpy.polynomial.legendre import leggauss

def normal(x,t):
    return np.exp(-x*x/(2*t))/math.sqrt(2*math.pi*t)

nodes,weights=leggauss(160)
times=[.12,.35,.7,1.]
radii=[.4,.68,.92,1.1]
x=nodes*radii[0]
w=weights*radii[0]
q=normal(x,times[0])
for i in range(1,len(times)):
    y=nodes*radii[i]
    q=normal(y[:,None]-x[None,:],times[i]-times[i-1])@(q*w)
    x=y
    w=weights*radii[i]
ratio=q/normal(x,1.)
pos=ratio[x>0]
b=radii[-1]
v=1-2*b*normal(b,1.)/(2*ndtr(b)-1)
half=[]
for z in (8.,16.,32.,64.,128.):
    # u=z(b−y) avoids cancellation in both exp(z*y) and the reflected kernel.
    def integ(u):
        s=u/z
        return math.exp(-u)*normal(b-s,1.)*(-math.expm1(-2*b*s))*z
    val,err=quad(integ,0,min(2*b*z,80),epsabs=1e-13,epsrel=1e-13)
    half.append(dict(z=z,scaled_killed_mgf=float(val),reported_error=float(err)))
print(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,grid_times=times,grid_radii=radii,legendre_order=160,absolute_symmetry_error=float(max(abs(ratio-ratio[::-1]))),largest_positive_grid_increment=float(np.max(np.diff(pos))),conditional_variance=float(np.dot(q*w,x*x)/sum(q*w)),truncated_free_variance=float(v),killed_reference_limit=float(2*b*normal(b,1.)),halfline_mgf=half),indent=2))
