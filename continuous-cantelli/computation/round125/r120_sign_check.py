"""High precision calibration of the explicit monotone nonexact step obstruction."""
import json
import platform
import mpmath as mp
import sympy as sp

mp.mp.dps=70
y=mp.mpf('0.1')
tail=lambda v: mp.erfc(v/mp.sqrt(2))/2
cdf=lambda v: tail(-v)
g=lambda v: mp.exp(-v*v/2)/mp.sqrt(2*mp.pi)
limit=y*y*cdf(1)+(2*y-1)*g(1)
rows=[]
for B in (7,12,24,60,160,400):
    def value(z):
        prefix=mp.exp(-z*y)*(cdf(B-z)-cdf(y-z))-mp.exp(z*y)*(tail(y+z)-tail(B+z))
        far=mp.exp(-z*y)*tail(B-z)-mp.exp(z*y)*tail(B+z)
        return (prefix+mp.exp(-z*z/2)*far)/z
    z=mp.mpf(B-1)
    scaled=mp.exp(z*y)*(z*mp.diff(value,z,2)-mp.diff(value,z))
    rows.append({'B':B,'scaled_signed_derivative_order_2':mp.nstr(scaled,24)})
assert mp.mpf(rows[-1]['scaled_signed_derivative_order_2'])<0
u=sp.symbols('u',positive=True)
rayleigh=sp.integrate(u**3*sp.exp(-u*u/2),(u,0,sp.oo))/sp.integrate(u*sp.exp(-u*u/2),(u,0,sp.oo))
assert rayleigh==2
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__,
    decimal_precision=mp.mp.dps,rayleigh_ratio=str(rayleigh),scaled_limit=mp.nstr(limit,24),cases=rows,seed=None),indent=2))
