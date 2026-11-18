"""Weighted quadrature for highly oscillatory superquadratic renewal gaps."""
import json,math
from pathlib import Path
from scipy.integrate import quad

w=2*math.pi

def finite(b,alpha,k):
    n=b**k
    a=n**(2-alpha)
    h=(b-2)/2
    g=0.; z=0.; phase=-1.
    for ell in range(1,71):
        T=h/b**ell
        m=quad(lambda x: math.exp(-a*x**alpha),0,T,epsabs=1e-14)[0]
        t=quad(lambda x:math.expm1(-a*x**alpha)*math.cos(w*(T-x)),0,T,epsabs=1e-24)[0]
        z+=2**ell*m
        g+=2**ell*phase*t
        phase*=math.cos(math.pi*(b-1)/b**ell)
    limit=math.ceil(2*(65/a)**(1/alpha))/2
    lead=0.
    for s in range(k):
        T=h*b**s
        u=min(T,limit)
        g+=2**(-s)*quad(lambda x:math.expm1(-a*x**alpha),0,u,weight='cos',wvar=w,epsabs=1e-13)[0]
        z+=2**(-s)*quad(lambda x:math.exp(-a*x**alpha),0,u,epsabs=1e-13)[0]
        lead+=2**(-s)*(-a*alpha*T**(alpha-1)*math.exp(-a*T**alpha)*(-1 if b%2 else 1)/w**2)
    return dict(base=b,alpha=alpha,k=k,coefficient=g/z,coefficient_scaled=(g/z)/a**(2/alpha),
                raw_numerator=g,boundary_prediction=lead,boundary_ratio=g/lead)

rs=[finite(b,p,k) for b in (3,4) for p in (2.5,3.,4.) for k in (2,4,6,8,10)]
for r in rs:
    if r['k']==10:
        assert (r['coefficient']>0)==bool(r['base']%2)
        print(r)
Path("supercritical_audit.json").write_text(json.dumps(rs,indent=2)+"\n")
