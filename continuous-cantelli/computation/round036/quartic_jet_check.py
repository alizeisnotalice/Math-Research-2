"""Exact checks of the Hermite identities and numerical checks of the gamma ratio."""
from math import comb, factorial
import json
import platform
import mpmath as mp
import sympy as sp
x,u=sp.symbols('x u', real=True)
for k in range(1,25):
    H=sp.hermite_prob(k,x)
    assert sp.expand(sp.diff(H,x)-x*H/2-
        (k*sp.hermite_prob(k-1,x)-sp.hermite_prob(k+1,x))/2)==0
    assert k*k*factorial(k-1)+factorial(k+1)==(2*k+1)*factorial(k)
# For constant source in the exact quartic model, the off-diagonal Nth density
# coefficient alone is B_N. This also checks the exponential base 4 sqrt(d).
def moment(power):
    return sp.factorial2(power-1) if power%2==0 else 0
for n in range(2,28,2):
    polynomial=sp.Poly(x**(2*n)*sp.hermite_prob(n,x),x)
    direct=sum(coeff*moment(power) for (power,),coeff in polynomial.terms())/factorial(n)
    assert direct==comb(2*n,n)*moment(n)
mp.mp.dps=75
out=[]
for n in (16,32,64,128,256,1024,4096):
    C=mp.mpf(comb(2*n,n))/4**n
    gaussian=mp.power(2,mp.mpf(n)/2)*mp.gamma(mp.mpf(n+1)/2)/mp.sqrt(mp.pi)
    numerator=mp.power(n,mp.mpf(5)/4)*mp.gamma(mp.mpf(n)-mp.mpf('3.5'))/mp.sqrt(mp.factorial(n))
    ratio=numerator/(C*gaussian)
    out.append({'N':n,'gamma_ratio':mp.nstr(ratio,15),
                'gamma_ratio_scaled_by_N_to_5_over_2':mp.nstr(ratio*mp.power(n,mp.mpf('2.5')),16)})
print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
     'mpmath_digits':mp.mp.dps,'exact_hermite_norm_checks':24,
     'exact_constant_source_moment_checks':13,'gamma_checks':out},indent=2))
