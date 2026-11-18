"""Independent coefficient algebra and finite double-log flux calibration."""

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


X_VALUES = ('1e2','1e4','1e8','1e12','1e16')
VARIANTS = ('matched','constant_plus_one','log_plus_one')
KNOTS = (-140,-80,-40,-20,-10,-5,-2,-1,0,1,2,5,10,20,40,80,140)


def exact_checks():
    y,z = sp.symbols('y z',positive=True)
    F = y**(-sp.Rational(1,2))*sp.exp(-1/(2*y))
    divided = {}
    for n in range(1,6):
        p = sp.Poly(sp.simplify(sp.diff(F,y,n).subs(y,z**-2)*sp.exp(z*z/2)/z),z)
        moments = sum(c*sp.factorial2(deg-1) for (deg,),c in p.terms())
        normalized = -sp.cancel(moments/sp.factorial(n))
        assert normalized == -sp.factorial2(2*n-1)
        divided[n] = str(normalized)
    A,B,b,gamma,logX = sp.symbols('A B b gamma logX')
    mean_r_tail = A*logX/2+A/4+B/2
    flux_residual = (5*gamma/4)*logX+5*b/4+gamma/8-sp.Rational(3,2)
    poly = sp.Poly(mean_r_tail-flux_residual,logX)
    solution = sp.solve(poly.all_coeffs(),(A,B))
    numeric = {k:sp.simplify(v.subs(gamma,sp.Rational(5,2))) for k,v in solution.items()}
    assert numeric[A]==sp.Rational(25,4)
    assert numeric[B]==5*b/2-sp.Rational(11,2)
    return divided, {str(k):str(v) for k,v in solution.items()}, {str(k):str(v) for k,v in numeric.items()}


def calibrate(x_label, variant, digits, method):
    with mp.workdps(digits):
        X = mp.mpf(x_label)
        logX = mp.log(X)
        gamma = mp.mpf(5)/2
        b = mp.log(4/mp.sqrt(mp.pi))
        A = mp.mpf(25)/4+(variant=='log_plus_one')
        B = (5*b-11)/2+(variant=='constant_plus_one')
        a = gamma*logX+b
        r = (A*logX+B)/X
        d = 1/(2*(X+a+r))
        C = mp.sqrt(mp.pi/2)

        def tail_change(z):
            if z==0:
                return mp.mpf(0)
            inverse = z/X
            logy = logX-mp.log(z)
            ay = gamma*logy+b
            ry = inverse*(A*logy+B)
            return mp.expm1(mp.log1p(inverse*(ay+ry))/2-ry)

        tail_cuts = (0,1,2,4,8,16,32,64,128,192)
        def integrate_log_z(func):
            return mp.quad(lambda t: mp.exp(-t)*func(mp.exp(-t)),tail_cuts,method=method)

        tail = (1+integrate_log_z(tail_change))/X

        def F(y):
            return mp.exp(-1/(2*y))/mp.sqrt(y)

        Fd = F(d)
        logd = mp.log(d)

        def pot_change(s):
            y = mp.exp(s)
            if s<logd:
                log_ratio = -mp.log1p(y/d)/2+y/(2*d*(d+y))
                return Fd*mp.expm1(log_ratio)-F(y)
            log_ratio = -mp.log1p(d/y)/2+d/(2*y*(y+d))
            raw = F(y)*mp.expm1(log_ratio)
            return raw-Fd if s<=0 else raw

        cuts = tuple(sorted({mp.mpf(t) for t in KNOTS} | {logd}))
        zero_index = cuts.index(mp.mpf(0))
        pot_left = mp.quad(pot_change,cuts[:zero_index+1],method=method)
        # A single endpoint value has no weight. Neither method evaluates the endpoint.
        pot_right = mp.quad(pot_change,cuts[zero_index:],method=method)
        potential = (pot_left+pot_right)/C
        theory = 1/X+(-a+mp.mpf(3)/2)/X**2
        V = integrate_log_z(lambda z: -mp.expm1(-(z/X)*(A*(logX-mp.log(z))+B)))/X
        return tail,potential,theory,V


def main():
    polys,general,specific = exact_checks()
    rows = []
    with mp.workdps(72):
        max_cross = mp.mpf(0)
        for label,variant in itertools.product(X_VALUES,VARIANTS):
            X = mp.mpf(label)
            primary = calibrate(label,variant,72,'gauss-legendre')
            check = calibrate(label,variant,50,'tanh-sinh')
            cross = X*X*max(abs(a-b) for a,b in zip(primary,check))
            max_cross = max(max_cross,cross)
            assert cross<mp.mpf('1e-24'),cross
            tail,pot,theory,V = primary
            normalized = X**2*(tail-pot)
            rows.append({'X':label,'variant':variant,
                         'normalized_flux_gap_X_squared':mp.nstr(normalized,38),
                         'normalized_flux_gap_X_squared_over_logX':mp.nstr(normalized/mp.log(X),38),
                         'potential_error_X_squared':mp.nstr(X**2*(pot-theory),27),
                         'tail_flux':mp.nstr(tail,40),'frozen_prepotential':mp.nstr(pot,40),
                         'integrated_residual_V':mp.nstr(V,40),
                         'cross_method_error_times_X_squared':mp.nstr(cross,12)})
            print(label,variant,mp.nstr(normalized,15),flush=True)
        out = {'scope':'Exact finite coefficient algebra and floating-point frozen prepotential/tail tests for explicit profiles.',
               'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,
               'methods_digits':{'gauss-legendre':72,'tanh-sinh':50},
               'prepotential_J_coefficients_divided_by_sqrt_pi_over_two':polys,
               'general_coefficients':general,'gamma_5_over_2_coefficients':specific,
               'b_numeric':mp.nstr(mp.log(4/mp.sqrt(mp.pi)),45),
               'max_cross_method_error_times_X_squared':mp.nstr(max_cross,12),'rows':rows}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('flux-results.json')
    path.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
