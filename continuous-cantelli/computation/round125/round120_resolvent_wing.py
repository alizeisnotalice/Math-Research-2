"""Rational inequalities for the outside-wing D02 certificate."""
import json
import platform
import sympy as s
import mpmath as mp
u=s.Symbol('u',real=True)
h=s.Symbol('h',positive=True)
E=lambda n:sum((-u*u/2)**j/s.factorial(j) for j in range(n+1))
S=lambda a,n:sum(a**j*u**(2*j+1)/s.factorial(2*j+1) for j in range(n+1))
upper=s.integrate(E(8)*(S(2,5)+128*u**13/s.factorial(13)),(u,0,2))
lower=s.integrate(E(7)*S(4,5),(u,0,2))
v=s.Rational(7,6)
Jupper=s.Rational(72,100)*(v*v+2)/(v**3+3*v)
assert upper<s.Rational(136,100) and lower>s.Rational(209,100)
assert Jupper<s.Rational(48,100)
margin=1-s.Rational(4,5)*upper+s.Rational(4,25)*lower-Jupper/2
assert margin>s.Rational(64,10000)
mp.mp.dps=50
C=lambda x:x*x-4*x/3+mp.mpf(8)/15-mp.mpf(8)/9*mp.exp(-3*x/2)
def jensen_margin(y):
    I=lambda n:mp.quad(lambda u:u**n*mp.exp(-u*u/2-y*u),[0,mp.inf])
    return (1-I(3)/(3*I(1)))**2-2*I(5)/(45*I(1))
cut=mp.findroot(jensen_margin,(mp.mpf('.39'),mp.mpf('.42')))
pstar=1-mp.mpf(3)/(2*mp.sqrt(5))
lam=(1+2*pstar)/3
a=-2*mp.log(pstar)/3
assert (1-a)**2-lam**2<0
print(json.dumps(dict(all_cut_cap=mp.nstr(lam+a,25),cutoff=mp.nstr(cut,25),python=platform.python_version(),sympy=s.__version__,mpmath=mp.__version__,digits=50,
                      upper_I1=str(upper),lower_I2=str(lower),upper_J=str(Jupper),strict_margin=str(margin),
                      threshold=mp.nstr(mp.findroot(C,(mp.mpf('.96'),1)),25)),indent=2))
