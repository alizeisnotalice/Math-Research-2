import json
import platform

import mpmath as mp
import sympy as sp

z,d=sp.symbols('z d', positive=True)
assert sp.integrate((d-z),(z,0,d))==d**2/2
assert sp.integrate(z*(d-z),(z,0,d))==d**3/6
assert sp.integrate(z*(d-z)**2,(z,0,d))==d**4/12
x,A,B,w=sp.symbols('x A B w', real=True)
P=(x-A)*(B-x)
assert sp.expand(P-w*sp.diff(P,x,2)/2-((x-A)*(B-x)+w))==0


def calc(digits):
    with mp.workdps(digits):
        t=mp.mpf('.8')
        B=mp.mpf('1.7')
        cdf=lambda x:mp.erfc(-x/mp.sqrt(2))/2
        p=lambda x:mp.exp(-x*x/(2*t))/mp.sqrt(2*mp.pi*t)
        D=lambda d:2*cdf(B/mp.sqrt(t))-cdf((B-d)/mp.sqrt(t))-cdf((B+d)/mp.sqrt(t))
        records=[]
        for width in (mp.mpf('.05'),mp.mpf('.3'),mp.mpf('1.1')):
            direct=mp.quad(lambda u:p(B-u)-p(B+u),[0,width])
            formula=D(width)
            assert abs(direct-formula)<mp.mpf('1e-40')
            cap=width*width/(mp.sqrt(2*mp.pi*mp.e)*t)
            assert formula<cap
            records.append({'width':mp.nstr(width,18),'reflected_exact_mass':mp.nstr(formula,18),
                            'quadratic_mass_cap':mp.nstr(cap,18)})
        return {'digits':digits,'t':str(t),'B':str(B),'records':records}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'exact_checks':['reflected polynomial budgets','backward heat quadratic transform'],
    'runs':[calc(50),calc(75)]},indent=2))
