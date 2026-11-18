"""Independent exact checks of pure central inverse coefficients, folds, and critical ratio."""
from math import comb
import mpmath as mp
import sympy as sp
import json,platform
u,w,b,e,a=sp.symbols('u w b e a')
rho=2/(3*sp.sqrt(3))
s=sp.sqrt(3)/3
# Evaluate derivative of the perturbed critical value in e: the motion of the
# critical point drops out because F'_0=0.
for S,root,expected in ((1,sp.I*s,a/9),(-1,s,-a/9),(-1,-s,-a/9)):
    F=u+S*u**3*(1+a*e*u)
    assert sp.simplify(sp.diff(F,u).subs({e:0,u:root}))==0
    assert sp.simplify(sp.diff(F,e).subs(u,root)-expected)==0
# When linear perturbation is absent, the cubic e phase at the positive imaginary fold:
h3=sp.symbols('h3',real=True)
phase3=sp.simplify(h3*(sp.I*s)**6/(sp.I*rho))
assert sp.simplify(phase3-sp.I*h3*sp.sqrt(3)/18)==0
for N in (4,8,12,16):
    central=sp.diff((1+b*u**2)**(-N-1),u,N).subs(u,0)/sp.factorial(N)
    assert sp.expand(central-comb(3*N//2,N//2)*b**(N//2))==0
mp.mp.dps=65
c=mp.mpf('0.8')
dH=(mp.root(32,3)/3-1)/(4*c)
def gamma_ratio(n,K):
    N=mp.mpf(n)
    A=1/(4*c)+dH
    gamma=3*N/4+2-mp.mpf(K)/4
    logB=(mp.loggamma(3*N/2+1)-mp.loggamma(N+1)-mp.loggamma(N/2+1)
          +N*mp.log(c)/4+N*mp.log(2)/4+mp.loggamma(N/4+mp.mpf('.5'))-mp.log(mp.pi)/2)
    lg=(mp.mpf('1.25')*mp.log(N)-mp.loggamma(N+1)/2-(N+3)*mp.log(c)/2
        +mp.loggamma(gamma)-gamma*mp.log(A)-mp.log(4))
    return mp.exp(lg-logB)
rows=[]
for N in (64,256,1024,4096):
    C=mp.mpf(comb(2*N,N))/4**N
    rho_m=2/(3*mp.sqrt(3))
    dar=C*rho_m**(-N)/comb(3*N//2,N//2)
    rows.append(dict(N=N,Darboux_to_Lagrange=mp.nstr(dar,16),
        K12_ratio=mp.nstr(gamma_ratio(N,12),16),K16_scaled_by_N=mp.nstr(N*gamma_ratio(N,16),16)))
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,
   digits=mp.mp.dps,c=str(c),dH=mp.nstr(dH,22),imaginary_p3=str(phase3),rows=rows),indent=2))
