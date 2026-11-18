import json
import sys

import mpmath as mp
import sympy as sp


def main():
    r, u, k = sp.symbols("r u k", positive=True)
    t = r**2 / (2 * u)
    transformed = sp.simplify(
        sp.exp(-r**2 / (2 * t)) / sp.sqrt(2 * sp.pi * t) * t**(-1-k) * (-sp.diff(t, u))
    )
    expected = 2**k / sp.sqrt(sp.pi) * r**(-1-2*k) * u**(k-sp.Rational(1, 2)) * sp.exp(-u)
    residual = sp.simplify(transformed / expected - 1)
    assert residual == 0

    mp.mp.dps = 60
    kappa = mp.mpf(1) / 4

    def gaussian(x, var):
        return mp.exp(-x**2 / (2 * var)) / mp.sqrt(2 * mp.pi * var)

    def c(power):
        return 2**power * mp.gamma(power + mp.mpf("0.5")) / mp.sqrt(mp.pi)

    def power_transform(b, r):
        p0 = gaussian(r, b) if b > 0 else 0

        def integrand(u):
            if not u:
                return 0
            return 4 * (p0 - gaussian(r, b + u**4)) * u**(-1 - 4 * kappa)

        return mp.quad(integrand, [0, mp.mpf("0.2"), 1, 2, mp.inf])

    distance = mp.mpf(1)
    f_zero = -c(kappa) * distance**(-1-2*kappa)
    derivative = -(1 + kappa) * c(kappa + 1) * distance**(-3-2*kappa)
    zero_error = abs(power_transform(mp.mpf(0), distance) - f_zero)
    assert zero_error < mp.mpf("1e-55")
    derivative_checks = []
    for power in (6, 10, 14, 18):
        age = mp.mpf(2) ** (-power)
        ratio = (power_transform(age, distance) - f_zero) / (age * derivative)
        derivative_checks.append([str(age), mp.nstr(ratio, 24)])

    centers = [mp.mpf(0)]
    for index in range(1, 7):
        shift = 2 * mp.mpf(3)**(-index)
        centers = [center + sign * shift for center in centers for sign in (-1, 1)]

    reflection_checks = []
    for power in (2, 5, 8, 11):
        age = mp.mpf(2)**(-power)
        deviation = mp.sqrt(age / (1 + age))
        weights = [gaussian(a, 1 + age) for a in centers]
        left = mp.fsum(
            weight * mp.erfc(a / ((1 + age) * deviation * mp.sqrt(2))) / 2
            for weight, a in zip(weights, centers)
        )
        ratio = left / mp.fsum(weights)
        assert abs(ratio - mp.mpf("0.5")) < mp.mpf("1e-55")
        reflection_checks.append([str(age), mp.nstr(ratio, 24)])

    print(json.dumps({
        "python": sys.version.split()[0],
        "sympy": sp.__version__,
        "mpmath": mp.__version__,
        "decimal_precision": mp.mp.dps,
        "symbolic_zero_age_residual": str(residual),
        "endpoint_quadrature_error": mp.nstr(zero_error, 10),
        "single_source_distance": str(distance),
        "single_source_kappa": str(kappa),
        "endpoint_kernel": mp.nstr(f_zero, 24),
        "endpoint_age_derivative": mp.nstr(derivative, 24),
        "age_difference_ratios": derivative_checks,
        "reflection_check_atoms": len(centers),
        "reflection_occupancy_ratios": reflection_checks,
    }, indent=2))


if __name__ == "__main__":
    main()
