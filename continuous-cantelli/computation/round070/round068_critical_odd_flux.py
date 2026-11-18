"""Symbolic coefficients and numerical checks of the critical heat identities."""
import json
import platform
import mpmath as mp
import sympy as sp

A,B,gamma,b = sp.symbols('A B gamma b')
coefficients = sp.solve([A/2-5*gamma/4,A/4+B/2-5*b/4-gamma/8+sp.Rational(3,2)],(A,B))
assert sp.simplify(coefficients[B]-(5*b/2-gamma-3))==0
d = sp.Symbol('d', positive=True)
H = d**(-sp.Rational(3,2))*sp.exp(-1/(2*d))
assert sp.simplify(sp.diff(sp.log(H),d)*2*d**2 - (1-3*d))==0
mass = sp.integrate(H/(2*d),(d,0,sp.oo))
assert sp.simplify(mass-sp.sqrt(sp.pi/2))==0

mp.mp.dps=48
p=lambda x: mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
shift=mp.mpf(1)
d0=mp.mpf('0.23')
def pdf(q, y):
    return p(y/mp.sqrt(q))/mp.sqrt(q)
def diff_integrand(z):
    q=1+d0*z*z
    return pdf(q,shift-z)+pdf(q,shift+z)-2*pdf(q,z)
frozen=mp.quad(diff_integrand, [0,1,5,mp.inf])
log_mass=mp.quad(lambda z: mp.log(z)*(p(z-shift)+p(z+shift)-2*p(z)), [0,1,5,mp.inf])
pred=-mp.exp(-1/(2*d0))/mp.sqrt(d0)*log_mass/mp.sqrt(mp.pi/2)
assert abs(frozen-pred)<mp.mpf('1e-42')
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,
                      digits=mp.mp.dps,A=str(coefficients[A]),B=str(coefficients[B]),odd_layer_mass=str(mass),
                      frozen_d=str(d0),shift=str(shift),direct=mp.nstr(frozen,24),logarithmic=mp.nstr(pred,24)),indent=2))
