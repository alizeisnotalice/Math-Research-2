"""Check uniform remote-derivative algebra and local periodic-minimum models."""

import json
import sys

import mpmath as mp
import sympy as sp


def symbolic():
    d, h, s, w = sp.symbols("d h s w", positive=True)
    ratio_exponent = sp.expand(((d-h)**2-d**2)/(2*s))
    assert sp.simplify(ratio_exponent + h*(d-h/2)/s) == 0
    derivative_factor = w*sp.exp(-h*w/2)
    stationary = sp.solve(sp.diff(derivative_factor, w), w)
    assert stationary == [2/h]
    maximum = sp.simplify(derivative_factor.subs(w, stationary[0]))
    assert maximum == 2*sp.exp(-1)/h
    return str(ratio_exponent), str(maximum)


def density(var, value):
    return mp.exp(-value**2/(2*var))/mp.sqrt(2*mp.pi*var)


def quadratic_profile(eps, radius):
    scale = mp.sqrt(2*eps)
    endpoint = mp.asinh(radius/scale)
    def integrand(t):
        x = scale*mp.sinh(t)
        variance = eps+2*mp.sin(x/2)**2
        return density(variance, x)*scale*mp.cosh(t)
    return 2*mp.quad(integrand, [0, endpoint/4, endpoint/2, 3*endpoint/4, endpoint])


def flat_profile(eps):
    scale = mp.sqrt(eps)
    def variance(x):
        return eps+(2*mp.sin(x/2)**2)**2
    core = mp.quad(lambda t: scale*density(variance(scale*t), scale*t), [-16, -4, 0, 4, 16])
    separated = mp.quad(lambda x: density(variance(x), x), [mp.mpf(2), mp.mpf(3)])
    return core, separated


def main():
    ratio_exp, max_expr = symbolic()
    with mp.workdps(100):
        rows = []
        for power in (20, 60, 120, 240):
            eps = mp.power(2, -power)
            q_value = quadratic_profile(eps, mp.mpf(1)/4)
            core, separated = flat_profile(eps)
            rows.append({
                "epsilon_power_two": -power,
                "quadratic_core_over_log": mp.nstr(q_value/mp.log(1/eps), 20),
                "flat_scaled_core": mp.nstr(core, 20),
                "flat_separated_positive_mass": mp.nstr(separated, 20),
                "flat_two_piece_lower_bound": mp.nstr(core+separated, 20),
            })
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "decimal_precision": mp.mp.dps,
            "remote_density_ratio_exponent": ratio_exp,
            "remote_derivative_maximum": max_expr,
            "quadratic_analytic_coefficient": mp.nstr(mp.exp(-1)/mp.sqrt(mp.pi), 20),
            "rows": rows,
        }, indent=2))


if __name__ == "__main__":
    main()
