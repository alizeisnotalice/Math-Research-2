import json
import platform

import mpmath as mp
import sympy as sp

u,a,q=sp.symbols('u a q',positive=True)
assert sp.simplify(2*sp.integrate(sp.exp(-a*u)*sp.cos(q*u),(u,0,sp.oo))-2*a/(a*a+q*q))==0
assert sp.integrate(2*a/(a*a+u*u),(u,-sp.oo,sp.oo))==2*sp.pi
r,d=sp.symbols('r d',positive=True)
assert sp.simplify(((1/r)-(1/(r+d)))*r/2-d/(2*(r+d)))==0


def compare(digits):
    with mp.workdps(digits):
        A=lambda u:1-3*u if u<=mp.mpf('.5') else u-1
        def transform(q):
            z=-q*q/2+1j*q
            if abs(z)<mp.mpf('1e-7'):
                return 2*mp.quad(lambda u:A(u)*mp.exp(-q*q*u/2)*mp.cos(q*u),[0,mp.mpf('.5'),1])
            def prim(u,m,b):
                return mp.exp(z*u)*((m*u+b)/z-m/z**2)
            return 2*mp.re(prim(mp.mpf('.5'),-3,1)-prim(0,-3,1)+prim(1,1,-1)-prim(mp.mpf('.5'),1,-1))
        vals=[]
        for eps in map(mp.mpf,('0.25','0.0625','0.00390625','0.0')):
            def physical(u):
                if u+eps==0:
                    return mp.inf
                return 2*A(u)*mp.exp(-u*u/(2*(u+eps)))/mp.sqrt(2*mp.pi*(u+eps))
            Ephysical=mp.quad(physical,[0,mp.mpf('.5'),1])
            Efourier=mp.quad(lambda q: mp.exp(-eps*q*q/2)*transform(q)/mp.pi,[0,1,5,20,mp.inf])
            assert abs(Ephysical-Efourier)<mp.mpf('1e-20')
            vals.append({'eps':str(eps),'energy':mp.nstr(Ephysical,18),'absolute_cross_error':mp.nstr(abs(Ephysical-Efourier),5)})
        return {'digits':digits,'values':vals}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'model':'finite continuous signed time density +/-1 on t<1/2,t>1/2 along x=1+t',
    'calibrations':[compare(50),compare(75)]},indent=2))
