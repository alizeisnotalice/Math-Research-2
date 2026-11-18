import json
import platform
from fractions import Fraction as F
from math import factorial
import mpmath as mp
mp.mp.dps=50
t=3*F(314,100)/8
partial=sum(t**j/factorial(j) for j in range(5))
assert partial>F(22,7)
g=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
def params(y):
    A=mp.quad(lambda u:u*g(y+u),[0,mp.inf])
    m=mp.quad(lambda u:u*u*g(y+u),[0,mp.inf])/A
    m2=mp.quad(lambda u:u**3*g(y+u),[0,mp.inf])/A
    d=mp.log(m/A)
    assert abs(m2-(2-y*m))<mp.mpf('1e-45')
    return m,m2,d
m,m2,d=params(1)
delta=3*mp.pi/8-mp.log(mp.pi)
kap=min(1+delta/d,2/m2)
for q in map(mp.mpf,['.1','1','3']):
    z=mp.sqrt(2*q)
    left=mp.quad(lambda u:mp.sinh(z*u)/z*g(u),[0,mp.inf])/g(0)
    right=mp.quad(lambda s:mp.exp(q*(1-s*s)),[0,1])
    assert abs(left-right)<mp.mpf('1e-47')
for n in range(1,17):
    a=(mp.power(2,1/mp.mpf(n))-1)/4
    assert 2-(1+2*a)**n>0
    for k in [1,mp.mpf("1.5"),12]:
        for u in [mp.mpf(".1"),2,5]:
            phi=lambda s:mp.sinh(mp.sqrt(2*s)*u)/mp.sqrt(2*s)
            D=mp.diff(phi,k)/phi(k)
            assert mp.diff(phi,k,n)/phi(k)<=D**n*(1+mp.mpf("1e-45"))
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,
    decimal_digits=mp.mp.dps, orders_checked=list(range(1,17)), rational_exponential_lower_bound=str(partial),
    delta=mp.nstr(delta,20),kappa=mp.nstr(kap,20),
    q_strip_numerator=mp.nstr(mp.log(kap),20)),indent=2))
