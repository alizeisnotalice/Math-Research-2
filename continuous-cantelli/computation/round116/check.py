"""Check the optimizing exponents symbolically and calibrate the score tails on f=1."""
import json
import math
import platform

import scipy
import sympy as sp
from scipy.special import ndtr


t, r, u, s = sp.symbols('t r u s', positive=True)
checks = {}
for sign in (1, -1):
    gap = r - sign*u
    exponent = sign*s*u - s**2*gap/(2*t)
    optimum = sign*t*u/gap
    checks[str(sign)] = {
        'derivative_residual': str(sp.simplify(sp.diff(exponent, s).subs(s, optimum))),
        'value_residual': str(sp.simplify(exponent.subs(s, optimum) - t*u**2/(2*gap))),
        'quadratic_residual': str(sp.simplify(sp.diff(exponent, s, 2) + gap/t)),
    }
    assert set(checks[str(sign)].values()) == {'0'}

rows = []
for tilt, cutoff, level in [(2, 3, 0.4), (10, 11, 0.7), (10, 8, 2.2), (40, 30, 4.0), (25, 28, 2.8)]:
    upper_score_cut = cutoff-tilt
    # f=1 gives W=N(0,1), X=t+W. Difference of survival functions avoids loss in the upper tail.
    positive_actual = max(0.0, ndtr(-level) - ndtr(-upper_score_cut))
    negative_actual = ndtr(min(upper_score_cut, -level))
    positive_bound = math.exp(-tilt*level**2/(2*(cutoff-level)))
    negative_bound = math.exp(-tilt*level**2/(2*(cutoff+level)))
    assert positive_actual <= positive_bound and negative_actual <= negative_bound
    rows.append(dict(t=tilt, R=cutoff, u=level, positive_actual=positive_actual,
                     positive_bound=positive_bound, negative_actual=negative_actual, negative_bound=negative_bound))
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,sympy=sp.__version__,
                      symbolic=checks,gaussian_cases=rows, random_seed=None), indent=2))
