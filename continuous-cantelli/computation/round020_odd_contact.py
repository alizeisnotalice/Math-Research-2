"""High-precision quadrature for the first derivative at a mixed contact.

The local deficit equals c*u on the left and d*u**3 on the right for
0 < u < r. Integrals here refer only to that fixed neighborhood.
"""

import json
import platform
import sys
from pathlib import Path

import mpmath as mp


def geometric_lines(start, stop):
    """Split each positive interval by successive factors of at most eight."""
    points = [start]
    if not start:
        points.append(min(mp.mpf(1), stop))
    while points[-1] < stop:
        points.append(min(8 * points[-1], stop))
    return points


def integrate_contact(x, c, d, n, r, method):
    eps = mp.mpf(10) ** -n
    h = mp.sqrt(eps)
    q0 = 1 / mp.sqrt(2 * mp.pi)

    def source(offset):
        return q0 * mp.exp(-(x + offset) ** 2 / 2)

    def first_kernel(u, variance):
        return q0 * u * mp.exp(-u**2 / (2 * variance)) / variance**mp.mpf('1.5')

    # The left kernel is negative, with u the positive physical displacement.
    # Using u = z**2 regularizes its integrable eps=0 endpoint singularity.
    def left(z):
        u = z**2
        return -2 * h * z * source(-u) * first_kernel(u, eps + c * u)

    left_cuts = [mp.mpf(0), *geometric_lines(h / mp.sqrt(c), mp.sqrt(r))]
    left_total = mp.quad(left, left_cuts, method=method)

    # In the core put u = h*t. The resulting integrand has a regular limit.
    def right_core(t):
        u = h * t
        return source(u) * first_kernel(t, 1 + d * h * t**3)

    switch = min(r, mp.root(eps / d, 3))
    right_total = mp.quad(right_core, geometric_lines(mp.mpf(0), switch / h), method=method)
    if switch < r:
        # In the outer segment put w = 1/u. This resolves the e^{-const/u}
        # decay on a fixed scale, even when switch is very close to zero.
        def right_outer(w):
            u = 1 / w
            return h * u**2 * source(u) * first_kernel(u, eps + d * u**3)

        right_total += mp.quad(right_outer, geometric_lines(1 / r, 1 / switch), method=method)
    return left_total, right_total, q0 * source(0)


def main():
    rows = []
    max_difference = mp.mpf(0)
    for x, c, d in ((1, 1, 1), (0, 2, 1)):
        for n in (4, 8, 12, 16, 24):
            with mp.workdps(40):
                baseline = integrate_contact(x, c, d, n, mp.mpf('0.08'), 'tanh-sinh')
            with mp.workdps(55):
                left, right, limit = integrate_contact(x, c, d, n, mp.mpf('0.08'), 'gauss-legendre')
                difference = max(abs(ref - check) for ref, check in zip((left, right, limit), baseline))
                max_difference = max(max_difference, difference)
                assert difference < mp.mpf('1e-35')
                total = left + right
                values = {
                    'scaled_left': left,
                    'scaled_right': right,
                    'scaled_sum': total,
                    'predicted_limit': limit,
                    'error_to_limit': total - limit,
                    'cross_method_difference': difference,
                }
                row = {'x0': x, 'c': c, 'd': d, 'epsilon': f'1e-{n}'}
                row.update({key: mp.nstr(value, 36) for key, value in values.items()})
                rows.append(row)
    result = {
        'scope': 'Uncertified numerical quadrature of local first derivative only',
        'python': sys.version,
        'mpmath': mp.__version__,
        'platform': platform.platform(),
        'radius': '0.08',
        'methods': {'tanh-sinh': 40, 'gauss-legendre': 55},
        'max_cross_method_difference': mp.nstr(max_difference, 12),
        'rows': rows,
    }
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('odd-contact-results.json')
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(f'{len(rows)} rows, max cross-method difference {result["max_cross_method_difference"]}')


if __name__ == '__main__':
    main()
