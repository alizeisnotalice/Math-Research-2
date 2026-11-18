"""Independent checks of the finite constraints and the wrong-sign call direction."""
import json
import platform
import mpmath as mp
import sympy as sp

z = sp.Symbol("z")
P = sp.Poly(-(1-z)**4*(z-sp.Rational(1,2))**2,z)
J = 3
basis = [(z-1)**k for k in range(1,4)]
matrix = sp.Matrix([[sp.diff(b,z,r).subs(z,1) for b in basis] for r in (1,2)]
                   +[[b.subs(z,sp.Rational(1,2)) for b in basis]])
assert matrix.det() != 0
assert P.eval(1)==0 and P.eval(sp.Rational(1,2))==0
mp.mp.dps=70
lam=mp.mpf(1)/16
C=mp.mpf(2)
coef = [mp.mpf(str(P.nth(k).p))/int(P.nth(k).q) for k in range(P.degree()+1)]
closed=mp.sqrt(2/mp.pi)*sum(c*mp.sqrt(C-j*lam) for j,c in enumerate(coef))
def integrand(x):
    if not x: return mp.mpf(0)
    u=x*x*lam/2
    # Factor out exp(6u) so the infinite quadrature grid stays well-conditioned.
    a=-mp.expm1(-u)
    b=1-mp.exp(-u)/2
    return mp.exp(-C*x*x/2+6*u)*a**4*b**2/(x*x)
fourier=2/mp.pi*mp.quad(integrand,[0,1,4,mp.inf])
assert closed>0 and abs(closed-fourier)<mp.mpf('1e-62')
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,
                      digits=mp.mp.dps,lambda_parameter=str(lam),node=mp.nstr(mp.sqrt(2*mp.log(2)/lam),24),
                      matched_moment_orders=[2,4,6],polynomial=str(P.as_expr()),jacobian_basis_determinant=str(matrix.det()),
                      wrong_frozen_absolute_moment_closed=mp.nstr(closed,32),fourier_integral=mp.nstr(fourier,32),
                      tail_limit=str(P.eval(0))),indent=2))
