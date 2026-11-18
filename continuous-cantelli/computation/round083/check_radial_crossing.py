"""Check Fourier energy, crossing moments and logarithmic-scale algebra."""

from fractions import Fraction
import json
import sys

import mpmath as mp
import sympy as sp


def symbolic_checks():
    t, lam = sp.symbols("t lam", positive=True)
    omega = sp.symbols("omega", real=True)
    primitive = sp.exp((-lam+sp.I*omega)*t)/(-lam+sp.I*omega)
    assert sp.simplify(sp.diff(primitive, t)-sp.exp((-lam+sp.I*omega)*t)) == 0
    transform = sp.simplify(-2*sp.re(primitive.subs(t, 0)))
    assert sp.simplify(transform - 2*lam/(lam**2+omega**2)) == 0

    tau, d, eps = sp.symbols("tau d eps", positive=True)
    shifted_heat_ratio_exponent = sp.factor(d**2/(2*tau) - d**2/(2*(tau+eps)))
    assert sp.simplify(shifted_heat_ratio_exponent - eps*d**2/(2*tau*(tau+eps))) == 0
    assert sp.simplify(sp.Rational(1, 2) - shifted_heat_ratio_exponent.subs(d**2, tau) - tau/(2*(eps+tau))) == 0

    delta = sp.symbols("delta", real=True)
    kappa = delta/16
    central_gap = sp.expand((1+delta)*(2-kappa)**2/4 - (1+kappa/2))
    assert sp.expand(central_gap - (delta**3-63*delta**2+928*delta)/1024) == 0
    assert Fraction(928-63, 1024) > Fraction(1, 2)
    assert Fraction(1, 8)-Fraction(1, 16) > 0

    return {
        "temporal_kernel_transform": str(transform),
        "mollified_heat_kernel_exponent": str(shifted_heat_ratio_exponent),
        "origin_central_gap_polynomial": str(central_gap),
        "energy_dyadic_power_margin": "1/16",
    }


def normal(x):
    return mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)


def check_trial():
    # A continuous bounded centered radial monotone trial matching the first four moments.
    cutoff = mp.mpf(3)

    def ramp(x):
        return min(mp.mpf(1), max(mp.mpf(0), x-cutoff))

    bounds = [mp.mpf(0), cutoff, cutoff+1, mp.inf]
    prob = 2*mp.quad(lambda x: ramp(x)*normal(x), bounds)
    second = 2*mp.quad(lambda x: ramp(x)**2*normal(x), bounds)
    covariance = 2*mp.quad(lambda x: (x*x-1)*ramp(x)*normal(x), bounds)
    variance = second-prob**2
    scale = 2*covariance/variance

    def f(x):
        return 1+scale*(ramp(x)-prob)

    fmin = f(0)
    assert fmin > 0
    first = 2*mp.quad(lambda x: f(x)*normal(x), bounds)
    discrepancy = 2*mp.quad(
        lambda x: (x*x-1)*(f(x)-1)*normal(x), bounds
    ) - mp.quad(lambda x: (f(x)-1)**2*normal(x), bounds)
    assert abs(first-1) < mp.mpf("1e-50")
    assert abs(discrepancy) < mp.mpf("1e-50")
    evaluations = []
    for t in (mp.mpf("0.5"), mp.mpf(1), mp.mpf(2), mp.mpf(4)):
        raw_sign = 2*mp.quad(lambda x: mp.cosh(t*x)*(f(x)-1)*normal(x), bounds)
        mgf = 2*mp.quad(
            lambda x: mp.cosh(t*x)*mp.exp(-t*t*f(x)/2)*normal(x), bounds
        )
        assert raw_sign > 0
        evaluations.append({
            "t": mp.nstr(t, 4),
            "exponential_sign_integral": mp.nstr(raw_sign, 16),
            "full_transform": mp.nstr(mgf, 16),
        })
    return {"trial_minimum": mp.nstr(fmin, 17), "evaluations": evaluations}


def moment_gap():
    out = []
    for a in map(mp.mpf, ("0.01", "0.1", "1", "3", "8")):
        survival = mp.erfc(a/mp.sqrt(2))/2
        excess_over_one = a*(normal(a)/survival-a)
        assert excess_over_one > 0
        out.append({"radius": mp.nstr(a, 4), "excess_over_one": mp.nstr(excess_over_one, 18)})
    return out


def asymmetric_tail():
    A, B = mp.mpf("0.7"), mp.mpf(2)
    polynomial = lambda x: (x+A)*(x-B)
    rows = []
    for side, radius, other in ((1, B, A), (-1, A, B)):
        survival = mp.erfc(radius/mp.sqrt(2))/2
        actual = mp.quad(lambda x: polynomial(side*x)*normal(x), [radius, mp.inf])/survival
        expected = 1+other*(normal(radius)/survival-radius)
        assert abs(actual-expected) < mp.mpf("1e-65")
        rows.append({"side": side, "conditional_quadratic_moment": mp.nstr(actual, 18)})
    return rows


def main():
    exact = symbolic_checks()
    with mp.workdps(80):
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "precision": mp.mp.dps,
            "symbolic": exact,
            "sign_and_low_moment_trial": check_trial(),
            "future_tail_second_moment_margin": moment_gap(),
            "asymmetric_tail_quadratic_check": asymmetric_tail(),
        }, indent=2))


if __name__ == "__main__":
    main()
