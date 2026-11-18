import json
import platform

import mpmath as mp
import sympy as sp

p, k, lam, B = sp.symbols('p k lam B', positive=True)
assert sp.simplify(lam**2*k**2/2*sp.integrate(sp.exp(lam**2*k**2*p/2),(p,0,B))-(sp.exp(lam**2*k**2*B/2)-1))==0
n, c = sp.symbols('n c', real=True)
x = sp.symbols('x', real=True)
assert sp.simplify(sp.integrate(sp.exp(sp.I*n*sp.sqrt(c)*x-x*x/2),(x,-sp.oo,sp.oo)).rewrite(sp.erf)-sp.sqrt(2*sp.pi)*sp.exp(-c*n*n/2))==0


def run(digits):
    with mp.workdps(digits):
        L=mp.mpf(2)*mp.pi
        v0=mp.mpf('0.7')
        vmax=mp.mpf(4)
        phi=lambda z:mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)
        lower=phi(L/(2*mp.sqrt(v0)))/mp.sqrt(vmax)
        def kernel(z,v):
            return sum(phi((z+j*L)/mp.sqrt(v))/mp.sqrt(v) for j in range(-12,13))
        records=[]
        for v in (v0,mp.mpf(2),vmax):
            row={'v':str(v),'min_sample_over_lower':mp.nstr(min(kernel(j*L/128,v) for j in range(-64,65))/lower,18)}
            for n in (0,1,3):
                freq=2*mp.pi*n/L
                coeff=mp.quad(lambda x:mp.cos(freq*x)*kernel(x,v),[-L/2,-L/4,0,L/4,L/2])
                expected=mp.exp(-freq*freq*v/2)
                assert abs(coeff-expected)<mp.mpf('1e-45')
                row[f'fourier_{n}']=mp.nstr(coeff,18)
            records.append(row)
        return {'digits':digits,'uniform_lower_selected':mp.nstr(lower,18),'wrapped_fourier':records}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'symbolic':'marked exponential antiderivative and Gaussian Fourier multiplier',
    'torus_period':'2*pi','wrap_range':[-12,12],
    'runs':[run(digits) for digits in (50,75)]},indent=2))
