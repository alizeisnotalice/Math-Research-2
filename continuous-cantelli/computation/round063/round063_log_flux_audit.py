"""Independent exact kernel derivatives and double-log-tail coefficient audit."""
import json
import platform
import sympy as sp
import mpmath as mp

z,d=sp.symbols('z d',positive=True)
r=z/sp.sqrt(1+d*z*z)
H=-r*sp.exp(-r*r/2)
K=z*(1+d*z*z)**(-sp.Rational(3,2))*(r*r-1)*sp.exp(-r*r/2)
assert sp.simplify(z*sp.diff(H,z)-K)==0
kernel=sp.simplify(sp.diff(H,d).subs(d,0)/z)
primitive_mass=sp.integrate(kernel,(z,0,sp.oo))
assert sp.simplify(primitive_mass+sp.sqrt(sp.pi/2))==0
mp.mp.dps=55
# lambda=log(t): choose a=1/(2d)=lambda+2.5log(lambda)+b.
# Multiplying the flux tail by Lbar=log(L) gives sqrt(2)e^-b asymptotically.
b=mp.log(4/mp.sqrt(mp.pi))
def flux_scale(B):
    q=mp.mpf(B)
    def integrand(u):
        # x=log(t)=q/z; with u=q/x, integrate on fixed [0,1].
        if not u: return mp.sqrt(2)*mp.exp(-b)
        x=q/u
        return mp.sqrt(2)*mp.exp(-b)*mp.sqrt(1+(mp.mpf('2.5')*mp.log(x)+b)/x)
    return mp.quad(integrand,[0,1])
res=[dict(log_L=B,scaled_tail=mp.nstr(flux_scale(B),22)) for B in (30,300,3000,30000)]
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,
  precision=mp.mp.dps,log_primitive_first_derivative_mass=str(primitive_mass),
  kernel_derivative_identity=True,selected_b=mp.nstr(b,22),
  predicted_scaled_tail=mp.nstr(mp.sqrt(2)*mp.exp(-b),22),rows=res),indent=2))
