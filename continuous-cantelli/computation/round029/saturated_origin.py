"""Comoving quadrature of the saturated logarithmic origin packet.

The model is f(u)=u²/[2(ell + log(ell)/2 + log(2)/2)], ell=log(1/u).
Integration is over a centered |z|≤16 window for fixed r=e^-4.
"""

import json
import platform
import sys
from pathlib import Path

import mpmath as mp


LOG_FREQUENCIES = ('1e3', '1e4', '1e6', '1e8', '1e12', '1e16', '1e20', '1e24')
CUTS = [-16, -10, -6, -3, 0, 3, 6, 10, 16]


def packet(log_frequency, digits, method):
    with mp.workdps(digits):
        L = mp.mpf(log_frequency)
        half = mp.mpf('0.5')
        d = half * mp.log(2)

        def D(ell):
            return ell + half * mp.log(ell) + d

        def root_equation(ell):
            return ell - L + mp.log(2 * D(ell))

        def root_derivative(ell):
            return 1 + (1 + half / ell) / D(ell)

        ell0 = mp.findroot(root_equation, L - mp.log(2 * L), solver='newton', df=root_derivative)
        D0 = D(ell0)
        y0 = 2 * D0
        width = 2 * mp.sqrt(D0)
        gram_scale = mp.exp(-2 * L)

        def kernel(z):
            displacement = width * z
            scaled_displacement = displacement / y0
            delta_ell = -mp.log1p(scaled_displacement)
            ell = ell0 + delta_ell
            delta_D = delta_ell + half * mp.log1p(delta_ell / ell0)
            y = y0 + displacement
            exponent = -(displacement - 2*delta_D)**2 / (4*(D0 + delta_D))
            gram = -gram_scale * y*y / 2
            return mp.sqrt(ell / D0) * mp.exp(exponent + gram) / ((1 + scaled_displacement) * mp.sqrt(mp.pi))

        assert 0 < y0 + CUTS[0]*width < y0 + CUTS[-1]*width < mp.exp(L-4)
        mass = mp.quad(kernel, CUTS, method=method)
        return mass, ell0, D0, y0, root_equation(ell0)


def main():
    main_method = 'tanh-sinh'
    verification = 'gauss-legendre'
    rows = []
    with mp.workdps(90):
        d = mp.log(2)/2
        constant = (1+d)/2
        max_cross_difference = mp.mpf(0)
        for label in LOG_FREQUENCIES:
            L = mp.mpf(label)
            primary, ell0, D0, y0, residual = packet(label, 70, main_method)
            secondary, *_ = packet(label, 90, verification)
            difference = abs(primary-secondary)
            max_cross_difference = max(max_cross_difference, difference)
            mass_deficit = 1 - secondary
            values = {'mass': secondary, 'deficit': mass_deficit,
                      'L_times_deficit_div_log_L': L * mass_deficit / mp.log(L),
                      'inverse_L_constant_estimate': L * mass_deficit - mp.log(L)/4,
                      'inverse_L_constant_error': L * mass_deficit - mp.log(L)/4 - constant,
                      'ell0': ell0, 'D0': D0, 'y0': y0, 'saddle_equation_residual': residual,
                      'cross_method_difference': difference}
            row = {'L': label}
            row.update({key: mp.nstr(value, 45) for key, value in values.items()})
            rows.append(row)
            print(f'L={label:s}  mass={mp.nstr(secondary, 29)}  ratio={mp.nstr(values["L_times_deficit_div_log_L"], 17)}'
                  f'  constant={mp.nstr(values["inverse_L_constant_estimate"], 19)}')
        assert max_cross_difference < mp.mpf('1e-67'), max_cross_difference
        result = {'scope': 'Floating-point quadrature of a comoving local window; analytic bounds handle omitted positive mass outside it',
                  'python': sys.version, 'mpmath': mp.__version__, 'platform': platform.platform(),
                  'log_frequency_definition': 'L=log(t)', 'radius': 'exp(-4)',
                  'center_equation': 'ell0=L-log(2*D(ell0)); y0=2*D(ell0)',
                  'coordinate': 'y=y0+2*z*sqrt(D0)', 'z_interval_splits': CUTS,
                  'methods_digits': {main_method: 70, verification: 90},
                  'predicted_limit_L_times_deficit_div_log_L': '0.25',
                  'predicted_inverse_L_constant': mp.nstr(constant, 45),
                  'max_cross_method_difference': mp.nstr(max_cross_difference, 20), 'rows': rows}
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('saturated-origin-results.json')
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(f'max cross-method difference {result["max_cross_method_difference"]}')


if __name__ == '__main__':
    main()
