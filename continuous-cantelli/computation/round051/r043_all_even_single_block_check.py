"""Independent exact coefficient check and high-order entropy check."""
from fractions import Fraction as Q
from math import factorial
import json,platform
import mpmath as mp

def coefficient(j):
    a=3*Q(1,j)-(2*Q(1,j-1) if j>=2 else 0)
    a-=Q(2,2**j*j)
    if j>=2:a+=Q(1,2**(j-1)*(j-1))
    a+=Q(2**j,j+1)-3*Q(2**(j-1),j)
    if j>=2:a+=Q(2**(j-1),j-1)
    return a
assert coefficient(2)==Q(1,12) and coefficient(3)==Q(1,24)
for j in range(4,201):
    formula=(Q(j-3)*(1-Q(2**(j-1),j+1))+Q(1,2**(j-1)))/(j*(j-1))
    assert coefficient(j)==formula and -Q(2**(j-1),j*(j-1))<=formula<0
assert Q(49,18)**3<Q(84375,4096)
mp.mp.dps=90
rows=[]
for mval in (3,4,5,10,100,1000,10000):
    m=mp.mpf(mval); d=m-1
    lam=d/(2*m-1)*(4*m/(2*m-1))**(m/d)
    x=1/m; nu=1/(m*d)
    f=nu*(2*m*mp.log(m)-2*d*mp.log(d)+mp.log(2*nu)-1)+(mp.mpf('.5')-d*nu)*(1+mp.log(lam/(1-2*d*nu)))
    excess=f-(1+mp.log(lam))/2
    A=-(3-2*x)*mp.log1p(-x)+(2-x)*mp.log1p(-x/2)-1-(1-x)*(1-2*x)/(2*x)*mp.log1p(-2*x)
    assert abs(excess-A/d)<mp.mpf('1e-76') and excess>0
    rows.append(dict(m=mval,excess=mp.nstr(excess,18),cubically_scaled=mp.nstr(excess*m**3,18)))
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,digits=90,
    symbolic_coefficients_verified_through=200,rows=rows),indent=2))
