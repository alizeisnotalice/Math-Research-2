import json
import platform

import mpmath as mp
import sympy as sp

z, alpha = sp.symbols('z alpha', positive=True)


def moment_polynomial(k):
    d = k % 2
    hermite = sp.Poly(sp.hermite_prob(k+2, z), z)
    result = 0
    for (degree,), coefficient in hermite.terms():
        ratio = sp.prod(alpha + d + 1 + 2*j for j in range((degree-d)//2))
        result += coefficient*ratio
    return sp.Poly(result, alpha)


def primitive_polynomial(k, power):
    remaining = sp.Poly(z**power * sp.hermite_prob(k+2, z), z)
    result = 0
    while not remaining.is_zero:
        degree = remaining.degree()
        coefficient = remaining.LC()
        assert degree >= 2 and degree % 2 == 0
        remaining -= coefficient*sp.Poly(sp.hermite_prob(degree, z), z)
        result -= coefficient*sp.hermite_prob(degree-1, z)
    assert sp.expand(sp.diff(result, z)-z*result) == sp.expand(
        z**power*sp.hermite_prob(k+2, z)
    )
    return result


for k in range(1, 19):
    d = k % 2
    n = (k-d)//2
    actual = moment_polynomial(k)
    assert actual == sp.Poly(sp.prod(alpha-j for j in range(d, k+1, 2)), alpha)
    for m, power in enumerate(range(d, k+1, 2)):
        expected = 2**n * sp.factorial(m) * (-1)**(n-m) * sp.factorial(n-m)
        assert actual.diff().eval(power) == expected


def check(power, k, beta, digits):
    with mp.workdps(digits):
        beta = mp.mpf(beta)
        polynomial = primitive_polynomial(k, power)
        h = sp.lambdify(z, polynomial, 'mpmath')
        cutoff = -mp.mpf(80)
        upper = mp.mpf(4)
        nodes = [cutoff, -8, -4, -2, 0, 1, 2, upper]

        def primitive(s):
            value = mp.exp(s)
            return h(value)*mp.exp(-value*value/2)

        d = k % 2
        normal_moment = mp.power(2, mp.mpf(power+d-1)/2)*mp.gamma(mp.mpf(power+d+1)/2)
        expected_area = -normal_moment * int(moment_polynomial(k).diff().eval(power))
        computed_area = mp.quad(primitive, nodes)
        assert abs(computed_area/expected_area-1) < mp.mpf('1e-50')

        def signal(t):
            return t**beta*(1+t**mp.mpf('-0.25')*mp.sin(mp.pi*t/7))

        def conv(log_scale):
            return mp.quad(lambda s: primitive(s)*signal(log_scale-s), nodes)

        values = []
        for scale in (20, 100, 500, 2000):
            at_scale = conv(scale)
            geometric_difference = conv(2*scale)-at_scale
            values.append({
                'L': scale,
                'H_convolution_scaled_ratio': mp.nstr(at_scale/(expected_area*mp.power(scale, beta)), 16),
                'geometric_difference_scaled_ratio': mp.nstr(
                    geometric_difference/(expected_area*(mp.power(2, beta)-1)*mp.power(scale, beta)), 16
                ),
            })
        return {
            'alpha': power, 'k': k, 'beta': str(beta), 'digits': digits,
            'finite_log_window': [str(cutoff), str(upper)],
            'primitive_gaussian_polynomial': str(polynomial),
            'primitive_area': mp.nstr(computed_area, 18),
            'runs': values,
        }


print(json.dumps({
    'python': platform.python_version(), 'sympy': sp.__version__, 'mpmath': mp.__version__,
    'exact_ledger': {
        'k_range': [1, 18], 'moment_product_verified': True,
        'all_simple_root_derivatives_verified': True,
        'polynomials': {str(k): str(moment_polynomial(k).as_expr().factor()) for k in [3, 4, 7, 8, 9]},
    },
    'calibrations': [check(power, k, beta, digits)
                     for power, k, beta in [(4, 8, '0.5'), (3, 7, '-0.5')]
                     for digits in (60, 85)],
}, indent=2))
