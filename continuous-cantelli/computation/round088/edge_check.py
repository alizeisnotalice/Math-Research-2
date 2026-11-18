import json
import platform

import mpmath as mp
import sympy as sp

s,x,k,j = sp.symbols('s x k j', real=True)
z=x-k*s
w=sp.exp(-(k*k+j*j)*s/2-k*z)*sp.sin(j*z)
assert sp.simplify(sp.diff(w,s)-sp.diff(w,x,2)/2)==0
u=sp.symbols('u',positive=True)
assert sp.integrate(u*sp.erfc(u/sp.sqrt(2))/2,(u,0,sp.oo))==sp.Rational(1,4)


def run(digits):
    with mp.workdps(digits):
        phi=lambda x,t:mp.exp(-x*x/(2*t))/mp.sqrt(2*mp.pi*t)
        T=mp.mpf('.75')
        b=mp.mpf('.9')
        lim=b*phi(b,T)/T
        records=[]
        for power in (6,12,24,40):
            h=mp.power(2,-power)
            root=mp.sqrt(h)
            rate=mp.quad(lambda u:(phi(b-root*u,T)-phi(b+root*u,T))*mp.erfc(u/mp.sqrt(2))/root,[0,1,4,10,mp.inf])
            jump=mp.erfc(root/mp.sqrt(2*T))
            records.append({'power':power,'reflection_rate_over_bound':mp.nstr(rate/lim,18),
                            'finite_horizon_jump_kernel_mass':mp.nstr(jump,18)})
        return {'digits':digits,'limiting_reflection_rate':mp.nstr(lim,18),'records':records}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'exact_checks':['translated sine solves heat equation','normal reflection integral = 1/4'],
    'runs':[run(digits) for digits in (50,75)]},indent=2))
