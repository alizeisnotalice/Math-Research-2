"""Heat-ball surrogate, positive-wing saddle quadratures on a normalized window."""
import json
import platform
import mpmath as mp
import sympy as sp

s,c=sp.symbols("s c",positive=True)
v=-sp.log(2*sp.pi*c**2*s)
b=sp.sqrt(s*v)
assert sp.simplify((c*b/(2*s))/sp.diff(b,s) - c*v/(v-1))==0
w,X=sp.symbols("w X",positive=True)
rho=1/sp.sqrt(w)
assert sp.simplify((1-rho**2)*sp.log(X/sp.sqrt(w))
                   +sp.log(w)/2+rho**2*sp.log(rho)-(1-1/w)*sp.log(X))==0

mp.mp.dps=45
def surrogate_ratio(ell, rho):
    # Exact profile implicit equation v-log(v)=2log(1/y).
    # Integrate ±10 Gaussian standard deviations in the positive saddle.
    V=2*ell-mp.log(2*ell)-2*mp.log(rho)
    def integrand(z):
        u=rho*V+mp.sqrt(V)*z
        if u<=0:
            return mp.mpf(0)
        A=2*(ell-mp.log(u))
        value=A+mp.log(A)
        for _ in range(5):
            value -= (value-mp.log(value)-A)/(1-1/value)
        return mp.exp(rho*u-u**2/(2*value)-rho**2*V/2)
    integ=mp.quad(integrand,[-10,-5,0,5,10])/mp.sqrt(2*mp.pi)
    return integ*mp.sqrt(V/(2*ell))

rows=[]
for w0 in [2,4]:
    rr=1/mp.sqrt(w0)
    for L in [64,512,4096]:
        q=surrogate_ratio(L,rr)
        rows.append(dict(w=w0,log_t=L,ratio=mp.nstr(q,18)))
    assert abs(q-1)<mp.mpf('.001')
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,
    mpmath=mp.__version__,quadrature_decimal_digits=mp.mp.dps,
    gaussian_window=[-10,10],symbolic_flux=True,symbolic_constant=True,
    surrogate_ratios=rows),indent=2))
