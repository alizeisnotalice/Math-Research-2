"""Numerical verification of the pair-tilt derivative budget and translated frame."""
import json
import math
import platform
import numpy as np
import scipy
from scipy.special import logsumexp
from numpy.polynomial.legendre import leggauss

# Finite nonlattice calibration; numerical evidence only.
j=np.arange(12,dtype=float)
x=j*0.91+0.17*np.sin(1.79*j)
width=0.10*np.exp(-0.12*j)
eps=.4/(1+j)
alpha=1.25
nodes,w=leggauss(48)
pts=np.concatenate((nodes-1,nodes+1))*.5
qw=np.tile(w*.5,2)
u=width[:,None]*pts
logprior=np.log(width[:,None]*qw)-.5*(x[:,None]+u)**2-.5*math.log(2*math.pi)
v=eps[:,None]+abs(pts)**alpha
def tilted(t):
    tmp=logprior-.5*t*t*v
    logA=logsumexp(tmp,axis=1)
    p=np.exp(2*logA-logsumexp(2*logA))
    q=np.exp(tmp-logA[:,None])
    mean=np.sum(q*v,axis=1)
    m=p@mean
    var=p@((mean-m)**2+np.sum(q*(v-mean[:,None])**2,axis=1)/2)
    phase=np.exp(1j*t*u)
    b=np.sqrt(p)*np.sum(q*phase,axis=1)
    bp=np.sqrt(p)*np.sum(q*phase*(1j*u-t*(v-m)),axis=1)
    return m,var,b,bp,p@np.sum(q*u*u,axis=1)
rows=[]
for t in (3.,7.,15.,30.):
    m,V,b,bp,U=tilted(t)
    h=t*2e-5
    dm=(tilted(t+h)[0]-tilted(t-h)[0])/(2*h)
    db=(tilted(t+h)[2]-tilted(t-h)[2])/(2*h)
    rows.append(dict(t=t,pair_mean=m,variance=V,norm=np.linalg.norm(b),
                     mean_derivative_residual=abs(dm+2*t*V),
                     vector_derivative_residual=np.linalg.norm(db-bp),
                     derivative_budget_ratio=np.linalg.norm(bp)**2/(2*U+4*t*t*V)))
T=15.
delta=min(np.diff(x))
gap=x[:,None]-x
gram=T*np.sinc(T*gap/(2*np.pi))**2
eig=np.linalg.eigvalsh(gram)
bound=4*math.pi**2/(3*T*delta**2)
assert T-bound<=eig[0]<=eig[-1]<=T+bound
print(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                      min_center_distance=delta,triangle_halfwidth=T,gram_extremes=[eig[0],eig[-1]],
                      gershgorin_radius=bound,rows=rows),indent=2))
