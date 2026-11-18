import json
import platform

import mpmath as mp
import sympy as sp

ell, d, r, e, t = sp.symbols('ell d r e t', positive=True)
objective = d*(r-d)/2
assert sp.solve(sp.diff(objective, d),d) == [r/2]
assert sp.simplify(objective.subs(d,r/2)-r*r/8) == 0
assert sp.diff(objective,d,2) == -1
level = 1-d/t
band = 1+d/t
assert sp.simplify(2/level-1-band) == d*(d+t)/(t*(-d+t))
rate = (1/level-1)*(1-e)**2*t**2/2
assert sp.simplify(rate-d*(1-e)**2*t/2-d**2*(e-1)**2*t/(2*(-d+t))) == 0
assert sp.simplify(sp.limit(rate*(band**2-1),t,sp.oo)-d*d*(e-1)**2) == 0


def run(digits):
    with mp.workdps(digits):
        records = []
        ee = mp.mpf('0.1')
        for rr in (mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)):
            dd=rr/2
            aa=rr/10
            mass = (mp.erf(aa/mp.sqrt(2))-mp.erf(aa/(2*mp.sqrt(2))))/2
            for tt in (40,160,640):
                hh=1-dd/tt
                ll=1+dd/tt
                exponent=(1/hh-1)*(1-ee)**2*tt*tt/2
                ratio=mp.exp(-exponent*(ll*ll-1))
                log_width=mp.log(2*mp.pi)/2-exponent-mp.log1p(-ratio)
                records.append({
                    'plateau_length': str(rr),
                    'chosen_deficit': str(dd),
                    't': tt,
                    'gaussian_test_mass': mp.nstr(mass, 18),
                    'log_width_bound': mp.nstr(log_width,18),
                    'log_width_bound_per_t':mp.nstr(log_width/tt,18),
                })
        return {'digits':digits, 'records':records}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'exact_checks':['one-sided modulus optimization','admissible width-band scale','leading exponent','geometric ratio limit'],
    'calibrations':[run(digits) for digits in (50,75)]},indent=2))
