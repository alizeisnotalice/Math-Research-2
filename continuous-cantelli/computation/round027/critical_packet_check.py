"""Numerical check of the real-action gap in the critical Fourier packet."""
import json
import math
import platform

import scipy
from scipy.optimize import minimize_scalar, brentq

width = 0.005
rows = []
for n in (4,6,8,10,16,32):
    k = n-2
    theta = math.pi/(2*(n-1))
    b, g = math.sin(theta), math.cos(theta)
    A = (n-1)*b/n
    q, p = k*b/n, k*g/n
    t = (k*b/2)**(1/n)
    Dcrit = t**(2*(n-1))/k

    def rate(r, sign, ratio):
        x = sign*r
        U = r**n/n
        return ratio*Dcrit/r**k + ((1+2*width*q)*U +
                                 width*((x-p)**2-q*q)/2)/(1+2*width*U)

    def min_budget(ratio):
        return min(minimize_scalar(rate,bounds=(0.7*t,1.3*t),
                                   args=(s,ratio),method='bounded',
                                   options=dict(xatol=1e-14)).fun for s in (-1,1))

    ratio = brentq(lambda h: min_budget(h)-A,0.97,1.0,xtol=1e-14)
    rows.append(dict(n=n,A=A,Dcrit=Dcrit,real_min_gap=min_budget(1.0)-A,
                     local_numeric_threshold_ratio=ratio,
                     small_width_gap_slope=(t-p)**2/2,
                     direct_first_order_gap=math.sin(2*theta)-t**(n+1)))
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,
                      fixed_packet_width=width,rows=rows),indent=2))
