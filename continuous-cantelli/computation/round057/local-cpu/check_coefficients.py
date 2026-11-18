"""Verify the leading and first corrections of the tested one-sided jets."""
import json
from pathlib import Path

import sympy as sp

z, c, b, d, a = sp.symbols("z c b d a")
rows = []
for n in (3, 4, 5, 6, 8, 10, 12):
    variance = c*z**n*(1+b*z+d*z*z)
    for sign in (-1, 1):
        source = sp.series(sp.exp(-sign*a*z-z*z/2), z, 0, 5).removeO()
        jet = sp.expand(sp.diff(source*variance, z, 2)/2
                        + sp.diff(source*variance**2, z, 4)/8)
        constant = sp.factor(jet.coeff(z, n-2))
        slope = sp.factor(jet.coeff(z, n-1))
        assert sp.expand(constant-c*n*(n-1)/2) == 0
        cubic = 45*c*c if n == 3 else 0
        assert sp.expand(slope-c*n*(n+1)*(b-sign*a)/2-cubic) == 0
        rows.append({"order": n, "sign": sign, "constant": str(constant), "slope": str(slope)})
Path("jump-coefficients.json").write_text(json.dumps({"sympy": sp.__version__, "rows": rows}, indent=2)+"\n")
