import json
import platform

import mpmath as mp
import sympy as sp

s, t, f, V, G = sp.symbols('s t f V G', positive=True)
w = (-V+sp.I*sp.sqrt(V)*G)/t
assert w.subs(G,0) == -V/t
assert sp.simplify((sp.Poly(w*w,G).nth(0)+sp.Poly(w*w,G).nth(2))-(V*V-V)/t**2) == 0
assert sp.expand(-s*s/2+(s*s-t*t)*f/2)==sp.expand(((f-1)*s*s-f*t*t)/2)


def run(digits):
    with mp.workdps(digits):
        k=mp.mpf(2)
        phi=lambda x: mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
        a,b=mp.mpf(2),mp.mpf(3)
        target=phi(b)/(1+1j*k)
        compact=[]
        for t in (20,80,320):
            limit_cal=mp.quad(lambda y: mp.exp(-(1+1j*k)*y)*phi(b-y/t),[0,1,4,16,t*(b-a)])/target
            compact.append({'t':t,'scaled_plateau_real_ratio':mp.nstr(mp.re(limit_cal),18),
                            'scaled_plateau_imag_ratio':mp.nstr(mp.im(limit_cal),18)})
        eps=mp.mpf('0.25')
        def B(x):
            h=2*(1+eps*mp.cos(k*x*x/2))/x**3
            return -mp.expm1(-mp.log1p(h)/2-x*x*h/(2*(1+h)))
        micro=[]
        for n in (500,8000,128000):
            for phase in (0,1):
                t=mp.sqrt((4*n+2*phase)*mp.pi/k)
                width=1/t**2
                normalized=t*mp.quad(B,[t,t+width])/width
                micro.append({'n':n,'phase_pi':phase,'origin':mp.nstr(t,18),
                              't_times_subcritical_mean_B':mp.nstr(normalized,18)})
        return {'digits':digits,'compact_plateau_calibration':compact,'synthetic_micro_calibration':micro}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'compact_test':'M=2, lambda=2, plateau [2,3], scaled endpoint integral',
    'micro_test':'f=1+2*x^(-3)*(1+.25*cos(x^2)), window width x^(-2)',
    'exact_checks':['conditional shift first and second moments','compact exponential rate'],
    'calibrations':[run(digits) for digits in (50,75)]},indent=2))
