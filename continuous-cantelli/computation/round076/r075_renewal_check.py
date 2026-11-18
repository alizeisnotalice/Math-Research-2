"""Direct summation over the gap genealogy for two unequal contractions."""
import json
import math
import platform
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import gamma, gammainc

alpha = 1.6
def tube_factor(t,r1,r2,d,levels):
    gap = 1-r1-r2
    sig = math.exp(-t)
    coeff = 2*sig*gamma(1/alpha)/alpha
    contributions=[coeff]
    for level in range(levels):
        k=np.arange(level+1)
        length=gap*r1**k*r2**(level-k)
        counts=np.array([float(math.comb(level,int(i))) for i in k])
        contributions.append(coeff*np.dot(counts,gammainc(1/alpha,(length/(2*sig))**alpha)))
    return math.exp((1-d)*t)*math.fsum(contributions)

rows=[]
for r1,r2 in ((1/3,1/9),(0.27,0.41)):
    d=brentq(lambda x:r1**x+r2**x-1,0,1,xtol=1e-14)
    for t in (10.,20.,30.,40.):
        F=tube_factor(t,r1,r2,d,145)
        G=tube_factor(t,r1,r2,d,175)
        shifted=[tube_factor(t+math.log(r),r1,r2,d,175) for r in (r1,r2)]
        rows.append(dict(r1=r1, r2=r2, d=d,t=t,F=G, depth_discrepancy=abs(F-G),
                         one_cell_ratios=[v/G for v in shifted],
                         renewal_residual=G-r1**d*shifted[0]-r2**d*shifted[1]))
assert max(x['depth_discrepancy'] for x in rows)<1e-8
assert max(abs(x['renewal_residual']) for x in rows)<1e-12
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,numpy=np.__version__,
                      alpha=alpha,depths=[145,175],dimension_xtol=1e-14,deterministic=True,rows=rows),indent=2))
