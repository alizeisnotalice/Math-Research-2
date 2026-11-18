"""Symbolic optimization and model Gaussian upper-bound calibration."""

import json
import sys

import mpmath as mp
import sympy as sp


def symbolic_checks():
    d, theta, w = sp.symbols("d theta w", positive=True)
    gain = theta*((d+w)**2-w*w/(1-theta))/2
    optimum = d*(1-theta)/theta
    assert sp.simplify(sp.diff(gain, w).subs(w, optimum)) == 0
    assert sp.simplify(gain.subs(w, optimum)-d*d/2) == 0
    assert sp.simplify(sp.diff(gain, w, 2)+theta**2/(1-theta)) == 0
    r = sp.symbols("r", positive=True)
    old_gap = 2*(d-r)
    new_gap = d-r
    difference = sp.simplify(old_gap**2/4 - new_gap**2/2)
    assert sp.expand(difference-(d-r)**2/2) == 0
    return {
        "optimal_observation_offset": str(optimum),
        "optimal_exponent_constant": str(gain.subs(w, optimum).simplify()),
        "finite_radius_coefficient": str(difference),
    }


def normal(variance, x):
    return mp.exp(-x*x/(2*variance))/mp.sqrt(2*mp.pi*variance)


def wall_model():
    b, d, r, t0 = map(mp.mpf, ("0.7", "1.3", "0.15", "1.5"))
    D = 2*(d-r)
    rows = []
    limiting_prefactor = 8*b*normal(t0, b)/(t0*D**2*mp.sqrt(4*mp.pi))
    for power in (4, 8, 12, 20):
        h = mp.power(2, -power)
        def integrand(v):
            source_diff = normal(t0,b-h*v)-normal(t0,b+h*v)
            return source_diff*mp.exp(-D*v/2-h*v*v/4)/(h*mp.sqrt(4*mp.pi))
        normalized = mp.quad(integrand, [0, 1, 4, 16, mp.inf])
        rows.append({
            "h_power_two": -power,
            "old_heat_divided_by_h_to_three_half_and_gaussian": mp.nstr(normalized, 20),
        })
    return {
        "normalized_model_rows": rows,
        "analytic_limit": mp.nstr(limiting_prefactor, 20),
    }


def loglog_model():
    A = mp.mpf("0.845")
    C = mp.mpf(2)
    rows = []
    for L in (50, 500, 5000):
        h = A/(L-2*mp.log(L)+C)
        ratio = h*h*mp.exp(L-A/h)
        rows.append({"log_inverse_radius": L, "width_ratio": mp.nstr(ratio, 20)})
    return {
        "normalized_model_rows": rows,
        "analytic_limit": mp.nstr(A*A*mp.exp(-C), 20),
    }


def main():
    exact = symbolic_checks()
    with mp.workdps(70):
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "precision": mp.mp.dps,
            "symbolic": exact,
            "reflected_wall_model": wall_model(),
            "critical_loglog_model": loglog_model(),
        }, indent=2))


if __name__ == "__main__":
    main()
