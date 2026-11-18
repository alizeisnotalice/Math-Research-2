"""Symbolic and high-precision checks for the potential and signed window bounds."""

import json
import sys

import mpmath as mp
import sympy as sp


def exact_checks():
    x, s = sp.symbols("x s", positive=True)
    gaussian = sp.exp(-x**2/(2*s))/sp.sqrt(2*sp.pi*s)
    assert sp.simplify(sp.diff(gaussian, x, 2) - 2*sp.diff(gaussian, s)) == 0

    h, lam = sp.symbols("h lam")
    inverse = ((1+4*lam*h)**sp.Rational(1, 2)-1)/(2*lam)
    inner_primitive = sp.Poly(sp.series(inverse**5/5, h, 0, 10).removeO(), h, lam)
    averaged = sp.Integer(0)
    for (power_h, power_lam), coeff in inner_primitive.terms():
        if power_lam % 2 == 0:
            k = power_lam // 2
            average_lam = sp.factorial2(2*k-1)/sp.Integer(k+1)
            averaged += coeff * average_lam * h**power_h
    assert sp.expand(averaged - (h**5/5+2*h**7+55*h**9)) == 0
    return str(averaged)


def cdf(z):
    return mp.erfc(-z/mp.sqrt(2))/2


def reference_window():
    A = 1/mp.sqrt(2*mp.pi)
    rows = []
    for power in (4, 8, 12, 20):
        rho = mp.power(2, -power)
        # Exact loss of q=1_[0,1] when smoothed with standard deviation rho.
        error = 2*rho*(A - mp.exp(-1/(2*rho**2))*A + (1/rho)*cdf(-1/rho))
        rows.append({
            "rho_power": power,
            "leakage_over_rho": mp.nstr(error/rho, 20),
        })
    return rows, mp.nstr(2*A, 20)


def signed_boundary():
    nodes, weights = mp.gauss_quadrature(24, "legendre")
    weighted_nodes = [((node+1)/2, weight/2) for node, weight in zip(nodes, weights)]
    c = 1/mp.sqrt(2*mp.pi)
    baseline_limit = mp.quad(lambda z: z*(1-z)*cdf(-z), [0, 1])*c
    rows = []
    for power in (4, 8, 16, 32):
        h = mp.power(2, -power)
        for sign in (1, -1):
            contributions = []
            for z, wz in weighted_nodes:
                x = h+h**2*z
                shape = z*(1-z)
                source = c * mp.exp(-x*x/2)
                for tau, wt in weighted_nodes:
                    scaled_standard_deviation = mp.sqrt((1+h*z)**4+sign*tau*h*shape)
                    p_window = (
                        cdf(-z/scaled_standard_deviation)
                        -cdf(-(1/h+z)/scaled_standard_deviation)
                    )
                    contributions.append(wz*wt*sign*shape*p_window*source)
            rows.append({
                "h_power": power,
                "residual_sign": sign,
                "signed_leakage_over_h_to_seven": mp.nstr(mp.fsum(contributions), 20),
            })
    return rows, mp.nstr(baseline_limit, 20)


def main():
    inverse = exact_checks()
    with mp.workdps(70):
        leakage, leakage_limit = reference_window()
        signed, signed_limit = signed_boundary()
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "decimal_precision": mp.mp.dps,
            "conditional_inverse_averaged_primitive": inverse,
            "fixed_source_leakage": leakage,
            "fixed_source_leakage_limit": leakage_limit,
            "signed_boundary_model": signed,
            "absolute_signed_leakage_limit": signed_limit,
        }, indent=2))


if __name__ == "__main__":
    main()
