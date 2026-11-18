"""Evaluate the kernel core for w(x)=x^2/log(1/x)^k at r=exp(-L)."""
import json
import platform
import scipy
from math import exp, log, pi, sqrt
from scipy.integrate import quad
from scipy.special import ndtr

def scaled_core(L, k):
    sigma = L**(-k/2)
    # x=r(1+sigma*t) so the frozen kernel has unit bandwidth in t.
    def integrand(t):
        ratio = 1+sigma*t
        if ratio <= 0:
            return 0.
        relative_sigma = ratio*(1-log(ratio)/L)**(-k/2)
        z = abs(t)/relative_sigma
        return 2*relative_sigma*(exp(-z*z/2)/sqrt(2*pi)-z*ndtr(-z))
    edges = sorted(set([-1/sigma, 0., 7/sigma] + [float(t) for t in (-12,-6,-3,-1,1,3,6,12,32,80) if -1/sigma<t<7/sigma]))
    vals = [quad(integrand,a,b,epsabs=1e-13,epsrel=1e-12) for a,b in zip(edges,edges[1:])]
    return {'k': k, 'L': L, 'ratio': sum(v[0] for v in vals), 'reported_quadrature_error':sum(v[1] for v in vals)}

print(json.dumps({'python':platform.python_version(), 'scipy':scipy.__version__, 'cases':[scaled_core(L,k) for k in (0.5,1.,2.) for L in (64.,1024.,16384.,262144.)]},indent=2))
