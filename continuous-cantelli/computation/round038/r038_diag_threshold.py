"""Exact interior root brackets and finite-order absolute diagonal weights."""
import json
import math
import platform
from fractions import Fraction
import scipy
from scipy.optimize import brentq
import numpy as np
from scipy.special import gammaln

def poly(x):
    return 1-x-(1-x/2)**4
lo,hi=Fraction(9126219746157,10**13),Fraction(9126219746159,10**13)
assert poly(lo)>0>poly(hi)
root=brentq(poly,float(lo),float(hi),xtol=1e-15)
limit=.5*math.sqrt(root)*(1-root/2)
rows=[]
for k in (64,256,1024,4096,16384,65536):
    i=np.arange(1,k//2,dtype=float)
    j=2*i
    log_beta=i*math.log(2)+gammaln(k/2+1)-gammaln(k/2-i+1)-gammaln(2*k+1)+gammaln(2*k-j+1)
    rate=np.exp(log_beta/j)*np.sqrt(j/math.e)
    arg=int(np.argmax(rate))
    rows.append(dict(k=k,q_at_max=j[arg]/k,max_value=rate[arg],reciprocal=1/rate[arg]))
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,q_root=root,
                      exact_bracket=[str(lo),str(hi)],limit_value=limit,
                      critical_factor=1/limit,rows=rows),indent=2))
