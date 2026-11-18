"""Checks for the endpoint fork saddle and the slow-tail critical coefficient."""
import json
import mpmath as mp
import sympy as sp
import platform
mp.mp.dps=55
R,c,K,p=map(mp.mpf,('1.25','1.6','0.25','1'))
rows=[]
for lam in (50,200,800):
    q=mp.mpf(lam)**2/2
    saddle=(K/(q*c))**(mp.mpf(1)/(2*p))
    hi=mp.mpf('.15')/saddle
    def weighted(z):
        if z==0:
            return 0
        u=saddle*z
        return (saddle*mp.exp(lam*u-q*c*u**p)*(-mp.expm1(-q*mp.exp(-K/u**p)))
                *mp.exp(-(R+u)**2/2)/mp.sqrt(2*mp.pi)
                *(1-mp.exp(-2*lam*(R+u)))*mp.exp(lam*mp.sqrt(2*c*K)))
    boundaries=sorted(set([mp.mpf(0),mp.mpf('.5'),mp.mpf('.8'),mp.mpf(1),mp.mpf('1.2'),mp.mpf(2),hi]))
    corr=mp.quad(weighted,boundaries)
    rows.append([lam,mp.nstr(R-mp.sqrt(2*c*K)+mp.log(corr)/lam,20)])
gamma=sp.Integer(3)
B=sp.Integer(2)
k=B**(1+gamma)*gamma**gamma/(2*(1+gamma)**(1+gamma))
crit=(1+gamma)**(1+gamma)*gamma**(-gamma)*k*2**(-gamma)
assert sp.simplify(crit-(B/2)**(1+gamma))==0
v,P,C=sp.symbols('v P C')
maxcoef=sp.expand((B*P*v/2-v**2/2-C/2).subs(v,B*P/2))
assert sp.simplify(maxcoef-(B**2*P**2/8-C/2))==0
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__,
                      dps=55,randomness='none',R=str(R),c=str(c),K=str(K),p=str(p),
                      split_limit=mp.nstr(R-mp.sqrt(2*c*K),25),split_rate_values=rows,
                      gamma3_B2_k_exact=str(k),second_order_verification=str(maxcoef)),indent=2))
