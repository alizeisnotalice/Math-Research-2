"""Numerically check the integrated heat leakage, including a rough oscillatory graph."""
import json, platform
import numpy as np
from scipy import __version__
from scipy.special import ndtr, roots_legendre
from scipy.integrate import quad
z, weight = roots_legendre(48)
def leakage(h, alpha, q):
    scale = h ** (q / 2)
    def edge(t, center):
        u = center + t * scale
        if abs(u) > 2*h or u == 0:
            return 0.
        w = abs(u)**q
        delta = abs(u)**alpha * (0.65 + 0.33*np.sin(1/abs(u)))
        s = w + (z+1)*delta/2
        # Probability discrepancy in absolute value, written using stable normal tails.
        if 0 < u < h:
            d = ndtr(-u/np.sqrt(s)) + ndtr((u-h)/np.sqrt(s))
        elif u >= h:
            d = ndtr((h-u)/np.sqrt(s)) - ndtr(-u/np.sqrt(s))
        else:
            d = ndtr(u/np.sqrt(s)) - ndtr((u-h)/np.sqrt(s))
        return delta * np.dot(weight,d)/(2*h**alpha)
    mid = h/(2*scale)
    total = 0.
    for center,lo,hi in [(0,-2*h/scale,mid),(h,-mid,h/scale)]:
        for a,b in [(max(lo,-36),0),(0,min(hi,36))]:
            total += quad(edge,a,b,args=(center,),epsabs=2e-8,limit=240)[0]
    return total
rows=[(h, leakage(h,5,4)) for h in [1/8,1/12,1/16,1/24,1/32]]
print(json.dumps(dict(python=platform.python_version(),scipy=__version__,alpha=5,q=4,scaled_leakage=rows),indent=2))
