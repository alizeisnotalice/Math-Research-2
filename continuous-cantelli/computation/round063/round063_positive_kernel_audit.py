"""Gaussian finite-difference and half-line coefficient checks."""
import json
import platform

import mpmath as mp
import sympy as sp

t = sp.Symbol("t", positive=True)
F = t**(-sp.Rational(1, 2)) * sp.exp(-1/(2*t))
ratio = sp.factor(sp.diff(F, t, 2)/F)
assert sp.simplify(ratio - (3*t*t-6*t+1)/(4*t**4)) == 0
mp.mp.dps = 60
root_two_pi = mp.sqrt(2 * mp.pi)
def gaussian_second(x):
    return (x*x-1) * mp.exp(-x*x/2) / root_two_pi
def odd_shift(s):
    return mp.quad(lambda z: z*z*(gaussian_second(s-z)-gaussian_second(s+z))/2, [0, 1, mp.inf])
shifts = []
for s in (mp.mpf("0.25"), mp.mpf(1), mp.mpf(3)):
    numeric = odd_shift(s)
    exact = mp.erf(s/mp.sqrt(2))
    assert abs(numeric-exact) < mp.mpf('1e-50')
    shifts.append(dict(shift=str(s), odd_kernel_mass=mp.nstr(numeric, 23)))
print(json.dumps(dict(python=platform.python_version(), sympy=sp.__version__,
                      mpmath=mp.__version__, digits=mp.mp.dps, curvature=str(ratio),
                      odd_kernel_mass_audits=shifts), indent=2))
