"""Check the Rodrigues transform, Laguerre reciprocal sums and critical beta bound."""
import json
import platform
from fractions import Fraction as F
from math import comb,factorial
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import hyp1f1
import sympy as s

a,x=s.symbols('a x')
for k in (2,4,6,8):
    explicit=sum(s.Rational(comb(k,j)*factorial(k),factorial(k+j))*(a*x)**j for j in range(k+1))
    diff=s.expand(s.diff(x**(2*k)*s.exp(a*x),x,k)*s.exp(-a*x))
    assert s.expand(diff-s.Rational(factorial(2*k),factorial(k))*x**k*explicit)==0
for k in range(2,257,2):
    p1=F(k,k+1)
    p2=F(comb(k,2),(k+1)*(k+2))
    assert p1*p1-2*p2==F(k*(2*k+1),(k+1)**2*(k+2))
beta=[]
for k in (16,64,256,1024):
    m=k//2
    ratio=F(1)
    peak=ratio
    for j in range(m-1):
        ratio*=F(8*k*(m-j),(2*k-2*j)*(2*k-2*j-1))
        peak=max(peak,ratio)
    beta.append({'k':k,'max_beta_over_h_power':float(peak)})
def response(k,t):
    root=np.sqrt(k)
    def kernel(z,tilted):
        val=root+z
        if val<=0:return 0.
        weight=np.exp(k*np.log1p(z/root)-(val*val-k)/2)
        if tilted:return weight*np.exp(t*(val-2*root))*hyp1f1(-k,k+1,-t*val)
        return weight
    lim=min(13.,root)
    den=quad(lambda z:kernel(z,False),-lim,13.,epsabs=1.e-12)[0]
    num=quad(lambda z:kernel(z,True),-lim,13.,epsabs=1.e-12)[0]
    return float(num/den)
rows=[{'k':k,'T(-1)':response(k,-1.),'T(+1)':response(k,1.)} for k in (32,128,512,1024)]
print(json.dumps({'python':platform.python_version(),'sympy':s.__version__,'scipy':scipy.__version__,
    'Rodrigues_polynomial_checks':4,'exact_reciprocal_sum_checks':128,'beta':beta,
    'quadrature':'chi-mode offset in [-min(sqrt(k),13),13]; illustration only','signed_response':rows},indent=2))
