import json
import platform

import mpmath as mp
import sympy as sp

alpha=sp.symbols('alpha',positive=True)
beta=2/alpha-1
assert sp.simplify(3-3*alpha/2-beta-(2-alpha)*(3*alpha-2)/(2*alpha))==0
checks=[]
for a in (sp.Rational(1,2),sp.Rational(3,4),sp.Rational(3,2)):
    theta=sp.Rational(9,10)
    n=int(sp.floor((2/a-1)/(2*(1-theta))))+1
    assert 2*n*(1-theta)>2/a-1
    checks.append({'alpha':str(a),'theta':str(theta),'even_derivative_order':2*n})


def run(digits):
    with mp.workdps(digits):
        c=mp.mpf(1)
        cutoff=mp.mpf(1)/4
        records=[]
        for a in map(mp.mpf,('1.2','1.5','1.8')):
            b=2/a-1
            A=(2*c)**(mp.mpf('.5')-1/a)*mp.gamma(-b/2)/(a*mp.sqrt(2*mp.pi*c))
            scale=[]
            for n in (24,64,128):
                y=mp.power(2,-n)
                q=y**b
                high=mp.log(cutoff/y**(2/a))
                low=-4*n*mp.log(2)/(1-a/2)
                mid=mp.log(1/(2*q))
                splits=sorted({p for p in (low,-40,-12,-5,0,5,12,40,mid,high) if low<=p<=high})
                def normalized(logz):
                    z=mp.exp(logz)
                    return z**(1-a/2)*mp.exp(-q*q*z**(2-a)/(2*c))*mp.expm1((2*q*z-1)/(2*c*z**a))
                val=mp.quad(normalized,splits)/mp.sqrt(2*mp.pi*c)
                scale.append({'n':n,'ratio_to_limiting_coefficient':mp.nstr(val/A,18)})
            records.append({'alpha':str(a),'limiting_coefficient':mp.nstr(A,18),'values':scale})
        return {'digits':digits,'records':records}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'numeric_model':'one forward gap wing, v=u^alpha for u in (0,1/4)',
    'exact_checks':checks,'numeric_runs':[run(digits) for digits in (45,70)]},indent=2))
