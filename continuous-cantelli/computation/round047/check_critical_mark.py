"""Symbolic and high-precision checks for the critical marked equation."""

import json
import sys

import mpmath as mp
import sympy as sp


def symbolic_checks():
    k, v, h = sp.symbols("k v h", positive=True)
    a = (1 + sp.I*k)**2
    assert sp.expand(a - 1 - 2*sp.I*k + k**2) == 0
    primitive = -sp.exp(-a*v/2)/a
    assert sp.simplify(sp.diff(primitive, v) - sp.exp(-a*v/2)/2) == 0

    r = sp.symbols("r", nonnegative=True)
    cumulative_gamma = sp.integrate(v*sp.exp(-v), (v, 0, r))
    assert sp.simplify(cumulative_gamma - (1 - (1 + r)*sp.exp(-r))) == 0

    multiplier = sp.simplify(a.subs(k, sp.Rational(1, 2))**(-1))
    z = sp.symbols("z", real=True)
    period = 4*sp.pi
    carrier = sp.exp(sp.I*z/2)
    shallow = -sp.expand_complex(multiplier*sp.conjugate(carrier)).as_real_imag()[0]
    deep = sp.cos(z/2)
    shallow_coefficient = sp.simplify(sp.integrate(shallow*carrier, (z, 0, period))/period)
    deep_coefficient = sp.simplify(sp.integrate(deep*carrier, (z, 0, period))/period)
    assert sp.simplify(shallow_coefficient + multiplier*deep_coefficient) == 0
    resultant_amplitude = sp.expand_complex(2*(shallow_coefficient+deep_coefficient))
    assert resultant_amplitude != 0
    return {
        "M_infinity_at_half": str(sp.expand_complex(multiplier)),
        "marginal_amplitude_at_half_for_unit_scale": str(resultant_amplitude),
        "gamma_two_cdf": str(cumulative_gamma),
    }


def quadrature_checks():
    multiplier_rows = []
    for depth in ("0.0001", "1", "10", "100"):
        for frequency in ("0.125", "0.5", "0.9", "0.99"):
            h, k = mp.mpf(depth), mp.mpf(frequency)
            a = (1+1j*k)**2
            d = -mp.expm1(-h/2)
            exact = -mp.expm1(-a*h/2)/(a*d)
            grid = [h*j/8 for j in range(9)]
            actual = mp.quad(lambda v: mp.exp(-a*v/2)/(2*d), grid)
            assert abs(actual-exact) < mp.mpf("1e-65")
            multiplier_rows.append({
                "H": depth,
                "k": frequency,
                "direct": mp.nstr(actual, 18),
                "absolute_error": mp.nstr(abs(actual-exact), 8),
            })

    shift = mp.mpf("0.3")
    lattice_rows = []
    for radius in (64, 256, 1024):
        signal = 2*mp.pi*mp.fsum(
            mp.sinc((2*mp.pi*n+shift)/6)**4 for n in range(-radius, radius+1)
        )
        lattice_rows.append({
            "sum_radius": radius,
            "signed_error": mp.nstr(signal-4*mp.pi, 18),
        })

    amplitude = mp.mpf("0.2")
    critical_length = mp.pi/2
    model_mass = mp.quad(lambda z: 1+amplitude*mp.cos(z), [0, critical_length])
    assert abs(model_mass - critical_length - amplitude) < mp.mpf("1e-70")
    return {
        "multiplier_quadratures": multiplier_rows,
        "unmarked_atom_lattice_test": lattice_rows,
        "width_pi_over_two_high_frequency_model_mass": mp.nstr(model_mass, 18),
        "width_pi_over_two_constant_mass": mp.nstr(critical_length, 18),
    }


def main():
    exact = symbolic_checks()
    with mp.workdps(80):
        numeric = quadrature_checks()
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "decimal_precision": mp.mp.dps,
            "symbolic": exact,
            "numeric": numeric,
        }, indent=2))


if __name__ == "__main__":
    main()
