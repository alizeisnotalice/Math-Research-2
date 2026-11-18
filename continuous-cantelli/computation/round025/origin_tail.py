"""High-precision calibration of compact-origin and tilted-tail asymptotics.

Run with Python, mpmath and SymPy. The optional argument is a JSON path.
The tail integrals use separately specified positive tail profiles.
"""

import json
import platform
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


R = '0.25'
CURVATURES = ('0.25', '1', '4')
LAMBDAS = (4, 16, 64, 256, 1024)
TAIL_LAMBDAS = (32, 128, 512, 2048)
TAIL_START = 8


def gaussian_cdf(x):
    return mp.erfc(-x / mp.sqrt(2)) / 2


def q0():
    return 1 / mp.sqrt(2 * mp.pi)


def model(c, family):
    if family == 'quadratic':
        return c, mp.mpf(0), mp.mpf(0)
    return c, -c/2, c/3


def intervals(lo, hi):
    candidates = (-32, -16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16, 32, 64)
    return [lo, *(v for v in candidates if lo < v < hi), hi]


def half_moments(c, tilt):
    """M_j=∫_0^∞ y^j exp(tilt*y-c*y²/2)dy for j=0,...,6."""
    first = mp.sqrt(2*mp.pi/c) * mp.exp(tilt**2/(2*c)) * gaussian_cdf(tilt/mp.sqrt(c))
    moments = [first, (1 + tilt*first)/c]
    for j in range(2, 7):
        moments.append(((j-1)*moments[-2] + tilt*moments[-1])/c)
    return moments


def core_coefficients(c, d, e, side, direction):
    moments = half_moments(c, side*direction)
    return tuple(q0()*term for term in (
        moments[0], -side*d*moments[3]/2,
        d*d*moments[6]/8 - e*moments[4]/2 - moments[2]/2,
    ))


def compact_half(c, d, e, frequency, side, direction, method):
    def integrand(t):
        x = side*t/frequency
        f_scaled = t*t*(c + d*x + e*x*x)
        return mp.exp(direction*side*t - f_scaled/2 - x*x/2)

    return q0() * mp.quad(integrand, intervals(0, mp.mpf(R)*frequency), method=method) / frequency


def quadratic_half(c, frequency, side, direction, compact):
    precision = 1 + c*frequency**2
    scale = mp.sqrt(precision)
    shift = side*direction*frequency/scale
    multiplier = mp.exp(mp.mpf(frequency**2)/(2*precision))/scale
    if not compact:
        return multiplier*gaussian_cdf(shift)
    r = mp.mpf(R)
    return multiplier*(gaussian_cdf(r*scale-shift) - gaussian_cdf(-shift))


def displayed(values):
    return {key: mp.nstr(value, 36) for key, value in values.items()}


def run_core(dps, method, frequencies):
    values = {}
    with mp.workdps(dps):
        for family in ('quadratic', 'cubic_quartic'):
            for c_label in CURVATURES:
                c, d, e = model(mp.mpf(c_label), family)
                for direction in (1, -1):
                    for frequency in frequencies:
                        values[family, c_label, direction, frequency] = tuple(
                            compact_half(c, d, e, frequency, side, direction, method) for side in (1, -1)
                        )
    return values


def verify_tail_symbolic():
    t, a, b, c, k, ell, m = sp.symbols('t A B C K L M')
    first_moment = k*t**3 + ell*t**4 + (m+6*k)*t**5
    second_moment = k*k*t**6 + 2*k*ell*t**7
    third_moment = k**3*t**9
    defect = first_moment/(2*t*t) - second_moment/(8*t**4) + third_moment/(48*t**6)
    target = a*t + b*t*t + c*t**3
    replacements = {k: 2*a, ell: a*a+2*b, m: 2*c-12*a+2*a*b+2*a**3/3}
    assert sp.expand(defect.subs(replacements)-target) == 0
    return [str(sp.expand(defect).coeff(t, j)) for j in range(1,4)]


def tail_deficit(params, frequency, method):
    k, ell, m = params

    def integrand(z):
        x = frequency + z
        residual = (k + ell/x + m/x**2)/x**3
        return q0()*mp.exp(-z*z/2)*(-mp.expm1(-frequency**2*residual/2))

    return mp.quad(integrand, intervals(TAIL_START-frequency, mp.inf), method=method)


def main():
    factors = verify_tail_symbolic()
    core = run_core(75, 'gauss-legendre', LAMBDAS)
    check = run_core(50, 'tanh-sinh', (LAMBDAS[0], LAMBDAS[-1]))
    with mp.workdps(75):
        max_cross = max(abs(x-y) for key in check for x,y in zip(core[key],check[key]))
        rows = []
        max_exact = mp.mpf(0)
        for (family, c_label, direction, frequency), integrals in core.items():
            c, d, e = model(mp.mpf(c_label), family)
            s = sum(integrals)
            for side, numeric in zip((1, -1), integrals):
                a,b,coeff3 = core_coefficients(c,d,e,side,direction)
                approximation = a/frequency + b/frequency**2 + coeff3/frequency**3
                values = {
                    'compact_integral': numeric, 'compact_half_fraction': numeric/s,
                    'limiting_half_fraction': gaussian_cdf(side*direction/mp.sqrt(c)),
                    'coefficient_lam_inv': a, 'coefficient_lam_inv_sq': b, 'coefficient_lam_inv_cube': coeff3,
                    'three_term_approximation': approximation,
                    'error_times_abs_lambda_fourth': mp.mpf(frequency)**4 * (numeric-approximation),
                }
                if family == 'quadratic':
                    exact = quadratic_half(c,frequency,side,direction,True)
                    max_exact = max(max_exact, abs(exact-numeric))
                    exact_unrestricted = quadratic_half(c,frequency,side,direction,False)
                    values.update({'closed_form_compact_integral': exact,
                                   'closed_form_unrestricted_halfline_integral': exact_unrestricted,
                                   'closed_form_unrestricted_halfline_fraction': gaussian_cdf(
                                       side*direction*frequency/mp.sqrt(1+c*frequency**2))})
                row = {'family': family, 'c': c_label, 'd': mp.nstr(d,20), 'e': mp.nstr(e,20),
                       'lambda': direction*frequency, 'side': side}
                row.update(displayed(values))
                rows.append(row)
        assert max_cross < mp.mpf('1e-40'), max_cross
        assert max_exact < mp.mpf('1e-55'), max_exact
        trows = []
        tail_max_cross = mp.mpf(0)
        for direction in (1,-1):
            c,d,e = model(mp.mpf(1), 'cubic_quartic')
            a,b,coeff3 = tuple(x+y for x,y in zip(core_coefficients(c,d,e,1,direction),
                                                  core_coefficients(c,d,e,-1,direction)))
            k = 2*a
            ell = a*a + 2*b
            m = 2*coeff3-12*a+2*a*b+2*a**3/3
            for frequency in TAIL_LAMBDAS:
                deficit = tail_deficit((k,ell,m), frequency, 'gauss-legendre')
                if frequency in (TAIL_LAMBDAS[0], TAIL_LAMBDAS[-1]):
                    with mp.workdps(50):
                        alternate = tail_deficit((k,ell,m), frequency, 'tanh-sinh')
                    tail_max_cross = max(tail_max_cross, abs(alternate-deficit))
                approx = a/frequency+b/frequency**2+coeff3/frequency**3
                center = sum(compact_half(c,d,e,frequency,side,direction,'gauss-legendre') for side in (1,-1))
                row = {'lambda': direction*frequency}
                row.update(displayed({'K':k,'L':ell,'M':m,'deficit_on_model_tail':deficit,
                                      'three_term_center_prediction':approx,'compact_origin':center,
                                      'deficit_minus_compact_times_abs_lambda_fourth':
                                          mp.mpf(frequency)**4*(deficit-center)}))
                trows.append(row)
        assert tail_max_cross < mp.mpf('1e-40'), tail_max_cross
        result = {'scope': 'Uncertified floating-point quadrature of explicit local and separate tail models, plus exact algebraic tail coefficient check',
                  'python': sys.version, 'mpmath': mp.__version__, 'sympy': sp.__version__,
                  'platform': platform.platform(), 'compact_radius': R, 'tail_start': TAIL_START,
                  'methods': {'gauss-legendre': 75, 'tanh-sinh_at_endpoint_frequencies': 50},
                  'tail_deficit_coefficients': factors,
                  'max_core_cross_method_difference': mp.nstr(max_cross,12),
                  'max_core_exact_quadratic_formula_error': mp.nstr(max_exact,12),
                  'max_tail_cross_method_difference': mp.nstr(tail_max_cross,12),
                  'compact_rows': rows, 'tail_rows': trows}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('origin-tail-results.json')
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(f'{len(rows)} compact rows and {len(trows)} tail rows; max core cross-method difference {result["max_core_cross_method_difference"]}')


if __name__ == '__main__':
    main()
