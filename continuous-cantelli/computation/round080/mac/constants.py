"""Exact rational bounds on the limiting Hermite tail constants."""
from fractions import Fraction as F
import json,math
from pathlib import Path
import mpmath as mp

def evaluation(poly,x):
    value=F(0)
    for p in reversed(poly):
        value=value*x+p
    return value

def square(poly):
    q=[F(0)]*(2*len(poly)-1)
    for i,v in enumerate(poly):
        for j,w in enumerate(poly): q[i+j]+=v*w
    return q

def integral(poly,x,power=0):
    return sum(v*x**(j+power+1)/(j+power+1) for j,v in enumerate(poly))

s=37
a,b=F('2.5128624172523393539654'),F('2.5128624172523393539655')
p=[F(0)]+[F(1,2**j*math.factorial(j+1)) for j in range(1,s+1)]
first=F(3,2)**(s+1)/math.factorial(s+2)
tail=first/(1-F(3,2)/(s+3))
assert evaluation(p,a)+tail<1<evaluation(p,b)
upper=p.copy();upper[0]+=tail
c_lo=a-integral(square(upper),a)
c_hi=b-integral(square(p),a)
der=[(j+1)*v for j,v in enumerate(p[1:])]
first_der=first*(s+1)/3
tail_der=first_der/(1-F(3,2)*(s+2)/((s+1)*(s+3)))
up_der=der.copy();up_der[0]+=tail_der
D_lo=2*(a-1)-2*integral(square(up_der),b,power=1)
D_hi=2*(b-1)-2*integral(square(der),a,power=1)
assert c_lo<c_hi and D_lo<D_hi
assert c_hi-c_lo<F(1,10**21) and D_hi-D_lo<F(1,10**20)
mp.mp.dps=62
x=mp.findroot(lambda x:mp.expm1(x/2)-x,2.5)
phi=lambda x:(mp.expm1(x/2)-x/2)/(x/2) if x else mp.mpf(0)
deriv=lambda x:mp.fsum(j*x**(j-1)/(2**j*mp.factorial(j+1)) for j in range(1,70))
C=mp.quad(lambda t:1-phi(t)**2,[0,x])
D=2*(x-1)-2*mp.quad(lambda t:t*deriv(t)**2,[0,x])
out=dict(root=[str(a),str(b)],certified_C_approx=[str(float(c_lo)),str(float(c_hi))],
         certified_D_approx=[str(float(D_lo)),str(float(D_hi))],C_bounds=[str(c_lo),str(c_hi)],
         D_bounds=[str(D_lo),str(D_hi)],mp_C=str(C),mp_D=str(D),mp_x=str(x))
Path('constants.json').write_text(json.dumps(out,indent=2)+"\n")
print(C,D)
