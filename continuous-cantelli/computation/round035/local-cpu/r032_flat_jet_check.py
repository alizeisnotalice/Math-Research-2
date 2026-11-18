"""Exact arithmetic checks for the jet estimate; slow-flat profile outside exp-flat classes."""
import json
import platform
from fractions import Fraction
import mpmath as mp
import sympy as sp

z,u,s=sp.symbols('z u s', real=True, positive=True)
# At u>z the heat potential is sqrt(s) phi((u-z)/sqrt(s))-(u-z)Phi((z-u)/sqrt(s)).
w=u-z
P=sp.exp(-w*w/(2*s))/sp.sqrt(2*sp.pi*s)
H=s*P-w*sp.erfc(w/sp.sqrt(2*s))/2
assert sp.simplify(sp.diff(H,z,2)-P)==0
for n in (2,3,4,5,6):
    assert sp.simplify(sp.diff(H,z,n)-sp.diff(P,z,n-2))==0
# K=1,m=2^-4k gives q≥m^-1/4/(2u) and q≥m^-1/2/2 by AM-GM.
# n≤m^-1/4/64 and log(4R/u)≤4/u then give n log≤q/8.
for k in range(6,51):
    N=2**(k-6)
    assert 4*N <= Fraction(2**k,16)
    assert 2**(2*k-1) == Fraction(2**(2*k),2)
# v=exp(-L log L) with L=log(1/u) is flat to all powers and violates every
# stretched-exponential upper envelope. q increases with L when log L>1.
mp.mp.dps=60
rows=[]
for L in [64,128,256,512]:
    ell=mp.mpf(L)
    log_inverse_m=ell*mp.log(ell)
    log_q=log_inverse_m-2*ell
    log_N=log_inverse_m/4-mp.log(64)
    rows.append({'log_inverse_R':L,'log_inverse_m':mp.nstr(log_inverse_m,12),
        'log_min_q':mp.nstr(log_q,12),'log_jet_cutoff_before_floor':mp.nstr(log_N,12)})
print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'symbolic_potential_orders_verified':[2,3,4,5,6],'dyadic_cutoff_cases':45,
    'profile':'exp(-L log L), L=log(1/u)','mpmath_digits':mp.mp.dps,'rows':rows},indent=2))
