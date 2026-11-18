"""Exact algebra for the antipodal budget and the two unequal real tilts."""
import json
import platform
import sympy as sp

D, k, h = sp.symbols('D k h', positive=True)
u = 1+h/sp.sqrt(D)
r = k*u**2
assert sp.simplify(D*(r-k)+2*k*D*(1-u)-k*h*h)==0
q=sp.Rational(3,2)*sp.log(3)-sp.Rational(5,2)
constant=sp.simplify(3+2*q)
assert sp.expand_log(constant,force=True)==3*sp.log(3)-2
eta=sp.sqrt((1-k)/(2*k))
g=sp.Symbol('g', positive=True)
logfactor=k*h*h+(1-k)*sp.log(g*(1-k)/(k*h))
assert sp.simplify(sp.diff(logfactor,h).subs(h,eta))==0
a,b,t=sp.symbols('a b t', positive=True)
S=(1+b*b)*t/(1-b)
a=(1-b)**2/(2*(1+b*b))
assert sp.simplify(a*S-t*(1-b)/2)==0
assert sp.simplify(a*S*S-t*t*(1+b*b)/2)==0
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,antipodal_constant=str(constant),
                      decimal=str(sp.N(constant,24)),two_tilt_exponent='k*eta**2',optimal_eta=str(eta),
                      power_likelihood_identity=True),indent=2))
