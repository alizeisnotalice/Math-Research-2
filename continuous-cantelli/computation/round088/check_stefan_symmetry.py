"""Check constitutive signs and the two odd-branch exponent scales."""

import json
import sys

import mpmath as mp
import sympy as sp


def exact_checks():
    c, a, b = sp.symbols("c a b", nonnegative=True)
    g = sp.symbols("g", positive=True)
    # Parametrize the only potentially mixed constitutive region.
    H1, H2 = g+a, g-b
    product_mixed = sp.factor((H1-H2)*a)
    assert product_mixed == a*(a+b)
    r, lam = sp.symbols("r lambda", positive=True)
    shift = sp.expand(-r*r/2+lam*r-lam**2/2)
    assert sp.expand(shift+(r-lam)**2/2) == 0
    plus = 3-8/(r+3)
    minus = 2-1/r
    sign_expr = sp.factor(plus-minus)
    assert sp.simplify(sign_expr-(r-1)*(r-3)/(r*(r+3))) == 0
    return {
        "mixed_enthalpy_product": str(product_mixed),
        "model_plus_minus_factorization": str(sign_expr),
        "model_limt_r_times_lower_tail_deficit": str(sp.limit(r*(2-minus),r,sp.oo)),
    }


def gaussian(x):
    return mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)


def synthetic_branches():
    def plus(x):
        return 3-8/(x+3)

    def minus(x):
        return 2-1/x

    rows = []
    # The branches agree on [0,1]. The particular common profile there
    # cancels from the odd identity, so it is not needed by the model.
    for lam in map(mp.mpf, (8,16,32,64,128)):
        def positive_integrand(y):
            x = 1+y/lam**2
            return (
                gaussian(x)*mp.exp(y/lam-lam**2*(plus(x)-1)/2)
                * (-mp.expm1(-2*lam*x))
                * (-mp.expm1(-lam**2*(minus(x)-plus(x))/2))
                / (2*lam**2)
            )

        def negative_integrand(z):
            x = lam+z
            return (
                gaussian(z)*mp.exp(-lam*z/(2*x))
                * (-mp.expm1(-2*lam*x))
                * (-mp.expm1(-lam**2*(plus(x)-minus(x))/2))/2
            )

        pos_end = 2*lam**2
        positive = mp.quad(positive_integrand, sorted(set([0,1,4,16,64,pos_end])))
        negative_bounds = sorted(set([3-lam, -min(lam/2,8), 0, 8, mp.inf]))
        negative = mp.quad(negative_integrand, negative_bounds)
        rows.append({
            "lambda": int(lam),
            "positive_log_exponential_rate_after_quadratic_removed": mp.nstr(
                1+mp.log(positive)/lam, 18
            ),
            "negative_log_exponential_rate_after_quadratic_removed": mp.nstr(
                mp.mpf("0.5")+mp.log(negative)/lam, 18
            ),
        })
    return rows


def main():
    exact = exact_checks()
    with mp.workdps(80):
        print(json.dumps({
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "decimal_precision": mp.mp.dps,
            "symbolic": exact,
            "synthetic_branch_rates": synthetic_branches(),
            "predicted_positive_rate": "1",
            "predicted_negative_rate": "0.5",
        }, indent=2))


if __name__ == "__main__":
    main()
