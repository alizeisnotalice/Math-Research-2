"""Independent finite LPs with exact Gaussian bin moments and direct quadrature."""
import json
import math
import warnings
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog, OptimizeWarning
from scipy.special import ndtr,ndtri
from scipy.stats import qmc

SQ=math.sqrt(2*math.pi)


def density(x):
    return np.exp(-.5*x*x)/SQ


def main():
    u,v,z,h=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(88789)).random_base2(5).T
    out=[]
    for x,w,rf,p in zip(u,v,z,h):
        l=-3+4*x
        r=l+.45+6*w
        m=(ndtr(r)-ndtr(l))*(.08+.83*rf)
        left=(.12+.76*p)*m
        a=ndtri(ndtr(l)+left)
        b=-ndtri(ndtr(-r)+(m-left))
        direct=[]
        for j in range(3):
            direct.append(quad(lambda x:x**j*density(x),l,a,epsabs=1e-14)[0]+
                          quad(lambda x:x**j*density(x),b,r,epsabs=1e-14)[0])
        mu=direct[1]
        peak=m+l*density(l)-a*density(a)+b*density(b)-r*density(r)
        assert abs(direct[0]-m)<1e-12 and abs(direct[2]-peak)<1e-12
        cert=-(a*b)*m+(a+b)*mu
        cert+=quad(lambda x:(x-a)*(x-b)*density(x),l,a,epsabs=1e-14)[0]
        cert+=quad(lambda x:(x-a)*(x-b)*density(x),b,r,epsabs=1e-14)[0]
        assert abs(cert-peak)<1e-12
        rows=[]
        for n in (64,128,256,512):
            e=np.linspace(l,r,n+1)
            p0=np.diff(ndtr(e))
            p1=-np.diff(density(e))
            p2=p0-np.diff(e*density(e))
            assert np.all(p0>0)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',category=OptimizeWarning)
                lp=linprog(-p2/p0,A_eq=np.vstack((np.ones(n),p1/p0)),b_eq=[1,mu/m],
                           bounds=np.column_stack((np.zeros(n),p0/m)),method='highs',
                           options={'threads':1})
            assert lp.success, lp.message
            estimate=-m*lp.fun
            assert estimate<=peak+2e-10
            rows.append(dict(bins=n,finite_lp=estimate,exact_gap=peak-estimate))
        assert np.all(np.diff([r['finite_lp'] for r in rows])>=-2e-10)
        out.append(dict(l=float(l),r=float(r),m=float(m),mu=float(mu),a=float(a),b=float(b),
                        exact_formula=float(peak),direct_quad=float(direct[2]),dual_quad=float(cert),lp=rows))
    Path('lp_audit.json').write_text(json.dumps(out,indent=2)+"\n")
    print("max quad",max(abs(x['direct_quad']-x['exact_formula']) for x in out))
    print("max dual",max(abs(x['dual_quad']-x['exact_formula']) for x in out))
    print("LP gaps by mesh", [max(x['lp'][i]['exact_gap'] for x in out) for i in range(4)])

if __name__=='__main__': main()
