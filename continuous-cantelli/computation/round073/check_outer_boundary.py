"""Check elementary heat identities and model boundary asymptotics."""

import json
import sys

import mpmath as mp
import sympy as sp


def symbolic():
    p, h = sp.symbols("p h", positive=True)
    exponents = [
        sp.factor((4/p-2)-(2/p-1)),
        sp.factor((3-3*p/2)-(2/p-1)),
        sp.factor(1-(2/p-1)),
    ]
    assert sp.simplify(exponents[0] - (2-p)/p) == 0
    assert sp.simplify(exponents[1] - (2-p)*(3*p-2)/(2*p)) == 0
    assert sp.simplify(exponents[2] - 2*(p-1)/p) == 0

    a, q, s = sp.symbols("a q s", positive=True)
    gaussian = sp.exp(-a*a/(2*s))/sp.sqrt(2*sp.pi*s)
    laplace = sp.integrate(sp.exp(-q*s)*gaussian, (s, 0, sp.oo))
    assert sp.simplify(laplace-sp.exp(-a*sp.sqrt(2*q))/sp.sqrt(2*q)) == 0
    return [str(x) for x in exponents]


def density(age, displacement):
    return mp.exp(-displacement**2/(2*age))/mp.sqrt(2*mp.pi*age)


def linear_model():
    rows = []
    c = mp.mpf("1.3")
    for y in map(mp.mpf, ("-1", "-0.1", "0", "0.1", "1")):
        endpoints = sorted(set([mp.mpf(0), abs(y)/10, abs(y), 1, 4, mp.inf]))
        val = mp.quad(lambda u: density(c*u, y+u), endpoints)
        predicted = mp.exp(-(y+abs(y))/c)
        assert abs(val-predicted) < mp.mpf("1e-35")
        rows.append({"y": mp.nstr(y, 5), "integral": mp.nstr(val, 20)})
    return rows


def comparable_cusp():
    rows = []
    p = mp.mpf("1.5")
    c = mp.mpf(1)
    rho = 1/p - mp.mpf("0.5")
    coefficient = (2*c)**(-1/p)*mp.gamma(1-rho)/(p*rho*mp.sqrt(mp.pi))
    delta = mp.mpf(1)/4
    for power in (8, 16, 24, 36):
        h = mp.power(2, -power)
        scale = h**(2/p)
        L = mp.log(delta/scale)

        def pieces(log_z):
            z = mp.exp(log_z)
            u = scale*z
            age = c*u**p
            baseline = density(age, u)
            main = -mp.expm1(-h*h/(2*age))*baseline
            cross = -baseline*mp.exp(-h*h/(2*age))*mp.expm1(u*h/age)
            return scale*z*main, scale*z*cross

        cut = sorted(set([mp.mpf(-60), mp.mpf(0), mp.log(h/scale), L]))
        main = mp.quad(lambda x: pieces(x)[0], cut)
        cross = mp.quad(lambda x: pieces(x)[1], cut)
        expected_scale = h**(2/p-1)
        rows.append({
            "h_power_two": -power,
            "difference_over_expected_scale": mp.nstr((main+cross)/expected_scale, 20),
            "cross_over_expected_scale": mp.nstr(cross/expected_scale, 20),
        })
    return rows, mp.nstr(coefficient, 20)


def logarithmic_saddle():
    out = []
    c = mp.mpf(1)
    for s in (32, 128, 512, 2048):
        center = s*mp.sqrt(c/2)
        lower = mp.log(4)
        nodes = [lower, center/2, center, 2*center, mp.inf]
        tail = mp.quad(lambda h: mp.exp(-(h-center)**2/h-s*mp.exp(-h)), nodes)
        log_integral = -2*center+mp.log(tail)
        out.append({"s": s, "negative_log_divided_by_s": mp.nstr(-log_integral/s, 20)})
    return out, mp.nstr(mp.sqrt(2), 20)


def main():
    exponents = symbolic()
    with mp.workdps(70):
        linear = linear_model()
        cusp, cusp_limit = comparable_cusp()
        saddle, saddle_limit = logarithmic_saddle()
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "decimal_precision": mp.mp.dps,
            "relative_error_exponents": exponents,
            "linear_models": linear,
            "subquadratic_superlinear_cusp": cusp,
            "analytic_cusp_coefficient": cusp_limit,
            "logarithmic_saddle": saddle,
            "analytic_logarithmic_saddle_rate": saddle_limit,
        }, indent=2))


if __name__ == "__main__":
    main()
