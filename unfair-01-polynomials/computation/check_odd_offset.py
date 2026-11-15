#!/usr/bin/env python3
"""Exact finite identity regression; the all-k theorem is proved analytically."""
import json
import sympy as sp

s = sp.symbols('s')
checks = 0
for k in range(1,50,2):
    W = sum((j+1)*(-s)**j for j in range(k+1))
    numerator = 1-(k+2)*s**(k+1)-(k+1)*s**(k+2)
    assert sp.expand((1+s)**2*W-numerator) == 0
    assert sp.expand((1-s**(k+1))*(1+s)-numerator
                     -s-(k+1)*s**(k+1)-k*s**(k+2)) == 0
    checks += 2
print(json.dumps({'exact_identity_checks':checks,'sympy':sp.__version__,
                  'scope':'odd k=1..49 regression only; analytic proof covers every odd k'}))
