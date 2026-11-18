"""Numerical calibration at c=1, x0=0, one critical negative block, g(0) normalized to 1."""
import json
import math
import platform
import scipy
from scipy.optimize import brentq

q = brentq(lambda t: 1-t-(1-t/2)**4,.6,.99)
type_bound=2/(math.sqrt(q)*(1-q/2))
def beta_ratio(i,k):
    return 2*(k//2-i)/((2*k-2*i)*(2*k-2*i-1))

def ratiosum(i,k,weighted):
    m=k//2
    h=1.
    out=1.
    for s in range(1,m-i):
        h *= beta_ratio(i+s-1,k)*(-.5)/s
        out+=h * ((m-i-s)/(m-i) if weighted else 1.)
        if abs(h)<1e-17:
            break
    return out

def coeff(k):
    m=k//2
    baseline=ratiosum(0,k,False)
    decrease=0.
    for l in range(2,k-1,2):
        if abs(l-q*k)>math.sqrt(k):
            continue
        i=l//2
        lb=i*math.log(2)+math.lgamma(m+1)-math.lgamma(m-i+1)-math.lgamma(2*k+1)+math.lgamma(2*k-l+1)
        term=math.exp(lb+l*math.log(type_bound)+math.lgamma(l+1)/2-1.75*math.log(l))
        decrease+=term*(m-i)*ratiosum(i,k,True)
    return dict(k=k,baseline=baseline,block_decrease_per_unit=decrease,
                diagonal_for_C_20=baseline-20*decrease)
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,critical_type=type_bound,
                      rows=[coeff(k) for k in (128,512,2048,8192,32768)]),indent=2))
