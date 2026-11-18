"""Independent exact power-series and fixed-order exponential-rate checks."""
import json
from fractions import Fraction as Q
import mpmath as mp
import sympy as sp
x=sp.symbols("x", positive=True)
ex=-(3-2*x)*sp.log(1-x)+(2-x)*sp.log(1-x/2)-1-(1-x)*(1-2*x)*sp.log(1-2*x)/(2*x)
series=sp.series(ex,x,0,15).removeO().expand()
assert series.coeff(x,2)==sp.Rational(1,12) and series.coeff(x,3)==sp.Rational(1,24)
for j in range(4,15):
    coeff=((j-3)*(1-Q(2**(j-1),j+1))+Q(1,2**(j-1)))/(j*(j-1))
    assert coeff==series.coeff(x,j) and 0 < -coeff < Q(2**(j-1),j*(j-1))
mp.mp.dps=75
rows=[]
for s in [3,4,5,6,8,16,32,64,128,256]:
    u=mp.mpf(1)/s
    # lambda_* = (s-1)/(2s-1) * (4s/(2s-1))^(s/(s-1))
    lam=(1-u)/(2-u)*(4/(2-u))**(1/(1-u))
    nu=u/(s-1)
    entropy=2*s*mp.log(2*s)-2*mp.log(2)-(2*s-2)*mp.log(2*s-2)
    r=1-2*(s-1)*nu
    gain=nu*(entropy+mp.log(2*nu)-1)+r*(1+mp.log(lam/r))/2-(1+mp.log(lam))/2
    A=-(3-2*u)*mp.log(1-u)+(2-u)*mp.log(1-u/2)-1-(1-u)*(1-2*u)*mp.log(1-2*u)/(2*u)
    assert abs(gain-A/(s-1))<mp.mpf('1e-68') and gain>0
    rows.append(dict(half_order=s,explicit_gain=mp.nstr(gain,24)))
print(json.dumps(dict(sympy=sp.__version__,mpmath=mp.__version__,rows=rows),indent=2))
