"""Cross-check the first real Laplace correction at critical equality."""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

from saddle import calc_constants, real_mp


def main():
    mp.mp.dps = 76
    n, z = sp.symbols("n z", positive=True)
    k = n-2
    J2 = (k+1+n-1)*z**(n-2)
    J3 = (-(k+1)*(k+2)+(n-1)*(n-2))*z**(n-3)
    J4 = ((k+1)*(k+2)*(k+3)+(n-1)*(n-2)*(n-3))*z**(n-4)
    lead = - (3*(n-1)**2+1)/(48*(n-1)*z**n)
    assert sp.simplify(5*J3**2/(24*J2**3)-J4/(8*J2**2)-lead) == 0
    assert 953**2 - 5*426**2 == 829
    assert sp.Rational(9,20)**8 < sp.Rational(13,20)**9
    rows = []
    for n in (4,6,8):
        k, theta, _, radius, D = calc_constants(n)
        H, b, x, c = map(mp.mpf, (".55","-.25",".9",".75"))
        zero = - (3*(n-1)**2+1)/(48*(n-1)*radius**n)
        for j in (4,8,14):
            s = mp.mpf(2)**-j
            lam = 2/(n*c*s**k)
            for sign in (-1,1):
                m, curv, _, _, val, tail = real_mp(n,k,D,radius,lam,s,sign*H,b,sign*x)
                amplitude = mp.exp(-sign*x*s*m-(s*m)**2/2)
                ratio = val/(mp.sqrt(2*mp.pi)*amplitude)
                rows.append({"n":n,"s_exponent":-j,"sign":sign,
                             "ratio_minus_one":mp.nstr(ratio-1,30),
                             "scaled_by_lambda":mp.nstr(lam*(ratio-1),25),
                             "expected_limit":mp.nstr(zero,25),
                             "scaled_residual":mp.nstr(lam*(ratio-1)-zero,20)})
    result = {"mpmath": mp.__version__, "dps": mp.mp.dps, "sympy": sp.__version__,
              "symbolic_first_correction": str(lead), "records": rows}
    Path("curvature-results.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
