"""Exact interval helpers extracted from the validated certificate generator."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
T=10**18
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def mul(a,b):
    c=[x*y for x in a for y in b];return min(c),max(c)
def scale(a,b):return mul(a,(b,b))
def inv(a):
    assert a[0]>0;return 1/a[1],1/a[0]
def arctan_inv(q,N):
    v=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),F(0))
    r=F(1,(2*N+1)*q**(2*N+1))
    return (v,v+r) if N%2==0 else (v-r,v)
# Machin identity: tan(2a)=5/12<1, hence 0<4a<pi/2;
# tan(4a-b)=1 with b=atan(1/239) and 0<4a-b<pi/2.
t2=2*F(1,5)/(1-F(1,5)**2);assert t2==F(5,12)<1
t4=2*t2/(1-t2*t2);assert t4==F(120,119)>F(1,239)
assert (t4-F(1,239))/(1+t4*F(1,239))==1
pi_machin=add(scale(arctan_inv(5,32),F(16)),scale(arctan_inv(239,8),F(-4)))
PI=(F(3141592653589793238,T),F(3141592653589793239,T))
assert PI[0]<pi_machin[0]<=pi_machin[1]<PI[1]
assert PI[1]<F(22,7)
@lru_cache(None)
def trig(a,kind):
    mid=(a[0]+a[1])/2
    # Rational centre rounded to a fixed denominator, with exact displacement.
    x=F((mid*T).__floor__(),T); rad=max(abs(a[0]-x),abs(a[1]-x))
    assert abs(x)<12
    term=x if kind=='sin' else F(1);value=term
    for j in range(1,41):
        denominator=(2*j)*(2*j+1) if kind=='sin' else (2*j-1)*(2*j)
        term*= -x*x/denominator;value+=term
    # Polynomial includes degree 81 (sin) or 80 (cos); degree-81 Taylor
    # polynomial for cos has its zero odd term. Derivatives are bounded by 1.
    remainder=abs(x)**82/F(factorial(82))+rad
    lo=value-remainder;hi=value+remainder
    # Outward rational rounding keeps downstream arithmetic small.
    return F((lo*T).__floor__(),T),F((hi*T).__ceil__(),T)
