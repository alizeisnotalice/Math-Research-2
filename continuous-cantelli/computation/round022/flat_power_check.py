"""Check Gaussian derivative asymptotics under a local flat power model."""

import json
import sys

import mpmath as mp
import sympy as sp


def moment_polynomial(k):
    alpha, x = sp.symbols("alpha x")
    parity = k % 2
    hermite = sp.Poly(sp.hermite_prob(k + 2, x), x)
    result = sp.Integer(0)
    for (power,), coefficient in hermite.terms():
        gaussian_moment_ratio = sp.prod(
            alpha + parity + 2 * j - 1 for j in range(1, (power - parity) // 2 + 1)
        )
        result += coefficient * gaussian_moment_ratio
    expected = sp.prod(alpha - root for root in range(parity, k + 1, 2))
    assert sp.expand(result - expected) == 0
    return str(sp.factor(result))


def hermite(k, z):
    previous, current = mp.mpf(1), z
    for j in range(1, k):
        previous, current = current, z * current - j * previous
    return previous if k == 0 else current


def scaled_kernel(k, s, t):
    z = s / mp.sqrt(t)
    return t**(-mp.mpf(k + 1) / 2) * hermite(k, z) * mp.exp(-z**2 / 2)


def normalized_scaled_difference(alpha, k, c_plus, c_minus, log_inverse_scale):
    # Normalize by phi(3)/sqrt(2pi), so the normalizing coefficient is unit-free.
    scale = mp.exp(-log_inverse_scale)
    h = scale ** (alpha - 2)
    cutoff = log_inverse_scale - mp.log(4)
    cut_points = sorted(set(range(-36, 37, 4)) | {cutoff})

    def integrand(log_s):
        s = mp.exp(log_s)
        u = scale * s
        core_u_alpha = min(mp.mpf(1) / 4, u) ** alpha
        baseline = scaled_kernel(k, s, 1)
        positive = scaled_kernel(k, s, 1 + c_plus * core_u_alpha / scale**2) - baseline
        negative = scaled_kernel(k, s, 1 + c_minus * core_u_alpha / scale**2) - baseline
        return (
            s * mp.exp(-u**2 / 2) / h
            * (mp.exp(-3 * u) * positive + (-1)**k * mp.exp(3 * u) * negative)
        )

    return mp.quad(integrand, cut_points)


def leading_coefficient(alpha, k, c_plus, c_minus):
    parity = k % 2
    gaussian_moment = 2 ** ((alpha + parity - 1) / 2) * mp.gamma((alpha + parity + 1) / 2)
    polynomial = mp.fprod(alpha - root for root in range(parity, k + 1, 2))
    return (c_plus + (-1)**k * c_minus) * gaussian_moment * polynomial / 2


def main():
    moments = {k: moment_polynomial(k) for k in range(11)}
    with mp.workdps(70):
        examples = [
            ("2.5", 1, 1, 0),
            ("2.5", 2, 1, 1),
            ("3", 2, 1, 1),
            ("4", 3, 1, 2),
            ("4", 4, 1, 2),
        ]
        results = []
        for alpha, k, c_plus, c_minus in examples:
            expected = leading_coefficient(mp.mpf(alpha), k, c_plus, c_minus)
            results.append(
                {
                    "alpha": alpha,
                    "derivative_order": k,
                    "c_plus": c_plus,
                    "c_minus": c_minus,
                    "leading_coefficient": mp.nstr(expected, 17),
                    "numerical": {
                        log_inverse_scale: mp.nstr(
                            normalized_scaled_difference(
                                mp.mpf(alpha), k, c_plus, c_minus, log_inverse_scale
                            ),
                            17,
                        )
                        for log_inverse_scale in (8, 16, 32)
                    },
                }
            )
        print(
            json.dumps(
                {
                    "python": sys.version.split()[0],
                    "mpmath": mp.__version__,
                    "sympy": sp.__version__,
                    "decimal_precision": mp.mp.dps,
                    "moment_polynomials": moments,
                    "results": results,
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
