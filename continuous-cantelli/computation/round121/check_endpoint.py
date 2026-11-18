import json
import sys

import mpmath as mp
import sympy as sp


def main():
    z, x, f = sp.symbols("z x f", real=True)
    tilt = sp.exp(z * x - z**2 * f / 2)
    ef_minus, ef_plus, cov, m, b, ey2, t = sp.symbols(
        "ef_minus ef_plus cov m b ey2 t", real=True
    )
    low_weighted = ef_minus + ef_plus + cov / m - 2 * (1 + t)
    low_after_normalization = low_weighted.subs(
        ef_plus, 2 * (1 + t + b * m - ey2) - ef_minus
    )
    symbolic = {
        "second_derivative": sp.simplify(
            sp.diff(tilt, z, 2).subs(z, 0) - (x**2 - f)
        ),
        "third_coefficient": sp.simplify(
            sp.diff(tilt, z, 3).subs(z, 0) / 6 - (x**3 / 6 - x * f / 2)
        ),
        "weighted_residual": sp.simplify(
            low_after_normalization - (2 * b * (m - b) + 2 * (b**2 - ey2) + cov / m)
        ),
    }
    assert all(value == 0 for value in symbolic.values())
    mp.mp.dps = 65

    def gaussian(y):
        return mp.exp(-y**2 / 2) / mp.sqrt(2 * mp.pi)

    def tail(y):
        return mp.erfc(y / mp.sqrt(2)) / 2

    def check(a, b):
        def integrand(y):
            return (y + a) * (y - b) * gaussian(y)

        tail_moment = mp.quad(integrand, [-mp.inf, -a]) + mp.quad(integrand, [b, mp.inf])
        excess = a * (gaussian(b) - b * tail(b)) + b * (gaussian(a) - a * tail(a))
        error = abs(tail_moment - tail(a) - tail(b) - excess)
        assert error < mp.mpf("1e-62")
        assert excess > 0
        return {"A": str(a), "B": str(b), "excess": mp.nstr(excess, 28), "error": mp.nstr(error, 5)}

    print(json.dumps({
        "python": sys.version.split()[0],
        "sympy": sp.__version__,
        "mpmath": mp.__version__,
        "decimal_precision": mp.mp.dps,
        "symbolic_residuals": {key: str(value) for key, value in symbolic.items()},
        "endpoint_cases": [check(mp.mpf(a), mp.mpf(b)) for a, b in (
            ("0.2", "0.2"), ("1", "1"), ("1.2", "0.8"), ("2", "3")
        )],
    }, indent=2))


if __name__ == "__main__":
    main()
