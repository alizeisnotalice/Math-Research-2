#!/usr/bin/env python3
"""Symbolic parameter identities for the analytic all-fixed-offset theorem."""
import json
import sympy as S

s,k,z,y = S.symbols('s k z y')
h=1+s
L=k+2+(k+1)*s
V=1+L*y
D=(6*k+3)*s**3+(k+4)*s**2-(4*k+9)*s+k+2
assert S.expand(s*(2*s-1)*V**2-y*h**2*((2*k+5)*s-k-2)
                -y*D-s*(2*s-1)*(L*y-1)**2)==0
assert S.factor(S.diff(D,s).subs(s,S.Rational(2,3))-(16*k+1)/3)==0
P=4*k**4-4*k**3-9*k**2-2*k+2
assert S.cancel(D.subs(s,k/(k+1))-P/(k+1)**3)==0
assert S.expand(P.subs(k,z+3))==4*z**4+44*z**3+171*z**2+268*z+131
# Independently check finite-sum cancellation; the proof covers all even k.
for n in range(4,42,2):
    W=lambda j:sum((i+1)*(-s)**i for i in range(j+1))
    assert S.expand((1-s)*W(n)+s*W(n+1)-W(n+2)
                     -s**(n+1)*((n+2)-(2*n+5)*s))==0
    assert S.expand(h*h*W(n)-1-(n+2)*s**(n+1)-(n+1)*s**(n+2))==0
print(json.dumps({'symbolic_parameter_identities':4,'finite_sum_regressions':38,
                  'sympy':S.__version__,
                  'scope':'all k identities plus even k=4..40 finite-sum regression; analytic proof supplies quantifiers'}))
