#!/usr/bin/env python3
"""Exact audit of the 7x18 structural escape and all 16 returned supports."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sympy as S
from audit_domination_failure import audit
from diagonal_domination import verify_certificate as verify_domination
from monomial_cycles import rows_for_support, find_certificate
from verify_monomial_cycles import verify as verify_cycle

I = [0, 2, 4, 7]
J = [0, 2, 4, 6, 8, 9, 10, 12, 13, 15, 18]
domination = audit(I, J)
variables, rows = rows_for_support(I, J)
w = [1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1]
assert len(rows) == 3 and len(w) == len(variables)
assert all(sum(a*b for a, b in zip(row['row'], w)) == 0 for row in rows)
assert all(n > 0 for n in w)

x = S.symbols('x')
a, b, c, d, e, f, g, h, i, j, k = S.symbols('a b c d e f g h i j k')
A = 1+a*x**2+b*x**4+x**7
B = 1+c*x**2+d*x**4+e*x**6+f*x**8+g*x**9+h*x**10+i*x**12+j*x**13+k*x**15+x**18
P = S.Poly(A*B, x)
equations = {q: S.expand(P.nth(q)-1) for q in [2, 9, 11, 4, 20, 13, 6]}
sub = {c: 1-a, g: a, d: 1-a*a, b: a*(2*a-1), j: 1-a,
       e: a*(1-a*(2*a-1))}
assert all(S.expand(equations[q].subs(sub)) == 0 for q in [2, 9, 11, 4, 20, 13])
p = S.expand(equations[6].subs(sub))
assert p == -5*a**3+4*a*a+a-1
alpha = (4+S.sqrt(31))/15
maximum = (-367+62*S.sqrt(31))/675
assert S.simplify(p.subs(a, alpha)-maximum) == 0
assert S.simplify(S.diff(p, a).subs(a, alpha)) == 0
assert 367**2-62**2*31 == 15525 > 0

here = Path(__file__).resolve().parent
batch = json.loads((here/'results/009-structural-search.json').read_text())
certs = []
for rec in batch['result']['records']:
    if 'I' not in rec:
        continue
    u, v = rec['I'], rec['J']
    assert u == sorted(set(u)) and v == sorted(set(v))
    assert u[0] == v[0] == 0 and u[-1] == 7 and v[-1] == 18
    assert len(u)+len(v) >= 11 and (set(u)&set(v))-{0}
    counts = Counter(ii+jj for ii in u for jj in v)
    assert sorted(t for t in counts if counts[t] == 1) == [0, 7, 18, 25]
    status = rec['verdict']['status']
    if status == 'diagonal_domination':
        verify_domination(u, v, rec['verdict']['certificate'])
        certs.append({'I': u, 'J': v, 'checked': 'direct_domination'})
    elif status == 'monomial_rule_not_certified_escape':
        cert = find_certificate(u, v)
        verify_cycle(cert)
        certs.append(cert)
    else:
        assert status == 'escapes_both_rules'
        assert (u, v) in [(I, J), (sorted(7-t for t in I), sorted(18-t for t in J))]
        certs.append({'I': u, 'J': v, 'checked': 'seven_diagonal_cubic_or_reversal'})
assert len(certs) == 16
print(json.dumps({'support_audit': domination, 'monomial_variables': variables,
                  'monomial_rows': rows, 'positive_right_kernel': w,
                  'equations': {q: str(eq) for q, eq in equations.items()},
                  'substitution': {str(k): str(v) for k, v in sub.items()},
                  'residual': str(p), 'maximum': str(maximum),
                  'all_16_supports_exactly_excluded': certs,
                  'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
