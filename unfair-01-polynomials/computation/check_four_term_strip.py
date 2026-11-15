#!/usr/bin/env python3
"""Exact symbolic identities and small odd-r Sturm certificate generation."""
import hashlib
import json
from pathlib import Path
import platform
import sympy as S
from verify_strip_sturm import verify

s, z, r, e = S.symbols('s z r e')
N = (1-z*z)*(1+s)**2-(1+e*s*z)*(1-e*(r+1)*z-e*r*s*z)
target = 2*s+s*s+e*((r+1)*z+(r-1)*s*z)-z*z+(r-1)*s*z*z+(r-1)*s*s*z*z
assert S.rem(S.expand(N-target), e*e-1, e) == 0
t = S.symbols('t')
br1 = 1-t*(1-s)
br2 = 1-s*br1-t*(1-s+s*s)
b2r1 = s-t*br1
R = (2*s-1)*(s-t)-t*t*(1-3*s+3*s*s)
assert S.expand(s*b2r1+t*br2+(1-s+s*s)-1-R) == 0
assert S.factor(R.subs(t, 1/(1+s))-s*(2*s**3+3*s*s-5*s+1)/(1+s)**2) == 0
a = S.Rational(5, 6)
K = 2+a-14*a**12-12*a**13-a**25
assert K > 0
assert all(x < 1 for x in [S.Rational(16, 14)*a*a, S.Rational(14, 12)*a*a, a**4])
assert (2*s**3+3*s*s-5*s+1).subs(s, a) == S.Rational(2, 27)

cases = []
for n in [3, 5, 7, 9, 11]:
    U = sum((-s)**i for i in range(n+1))
    b = [S.Integer(1)]
    for k in range(1, 2*n):
        b.append(S.expand(1-s*b[-1]-(U*b[k-n] if k >= n else 0)))
    H = S.Poly(S.cancel(b[-1]/(s*(s-1))), s)
    RR = S.Poly(R.subs(t, U), s)
    intervals = [(lo, hi) for (lo, hi), mult in H.intervals(eps=S.Rational(1, 10000))
                 if 0 < lo < hi < 1]
    assert len(intervals) == 1 and H.count_roots(0, 1) == 1
    lo, hi = intervals[0]
    c = {'r': n, 'H': [int(H.nth(i)) for i in range(H.degree()+1)],
         'R': [int(RR.nth(i)) for i in range(RR.degree()+1)],
         'interval': [str(lo), str(hi)]}
    c['independent_sturm'] = verify(c)
    cases.append(c)

here = Path(__file__).resolve().parent
print(json.dumps({'python': platform.python_version(), 'sympy': S.__version__,
                  'general_parity_identity': True, 'excess_identity': True,
                  'K13_at_5_6': str(K), 'small_odd': cases,
                  'sha256': {f: hashlib.sha256((here/f).read_bytes()).hexdigest()
                             for f in ['check_four_term_strip.py', 'verify_strip_sturm.py']}}, indent=2))
