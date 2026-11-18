import json
import sys

import mpmath as mp
import sympy as sp


def main():
    h, s, z = sp.symbols("h s z", positive=True)
    p = sp.exp(-z**2 / (2 * s)) / sp.sqrt(2 * sp.pi * s)
    next_h = h / (1 + h)
    delta = h + h**2
    symbolic = {
        "heat_age_derivative": sp.simplify(sp.diff(p, s) / p - (z**2 - s) / (2 * s**2)),
        "inverse_age_shift": sp.simplify(1 / h - 1 / delta - 1 / (1 + h)),
        "largest_age_in_bin": sp.simplify((delta - next_h) / h**2 - 1 - 1 / (1 + h)),
    }
    assert all(residual == 0 for residual in symbolic.values())

    mp.mp.dps = 70
    distance = mp.mpf(1)
    action = distance**2 / 2
    gamma_half = mp.sqrt(2 * mp.pi)

    def bin_fraction(h_value):
        # M(h_next) / M(h) evaluated without the common exp(-action/h).
        return 1 - (1 + h_value) ** (-mp.mpf(5) / 2) * mp.exp(-action)

    def wall_ratio(h_value):
        # Substitute b-x=h*u in the half-line source density b-x.
        return mp.quad(
            lambda u: u * mp.exp(-distance * u - h_value * u**2 / 2),
            [0, 1, 8, mp.inf],
        ) / gamma_half

    def record(h_value):
        return {
            "h": str(h_value),
            "normalized_bin_mass": mp.nstr(bin_fraction(h_value), 22),
            "normalized_wall_heat": mp.nstr(wall_ratio(h_value), 22),
        }

    data = {
        "python": sys.version.split()[0],
        "sympy": sp.__version__,
        "mpmath": mp.__version__,
        "decimal_precision": mp.mp.dps,
        "symbolic_residuals": {name: str(value) for name, value in symbolic.items()},
        "distance": str(distance),
        "records": [record(mp.mpf(2)**-k) for k in (8, 16, 24)],
        "normalized_bin_limit": mp.nstr(1 - mp.exp(-action), 22),
        "normalized_wall_limit": mp.nstr(1 / (gamma_half * distance**2), 22),
    }
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
