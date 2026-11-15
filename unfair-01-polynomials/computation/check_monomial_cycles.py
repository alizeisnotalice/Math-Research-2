#!/usr/bin/env python3
"""Reproduce four exclusions, two controls and adversarial certificate checks."""
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import scipy
from monomial_cycles import find_certificate, rows_for_support
from verify_monomial_cycles import verify

HERE = Path(__file__).resolve().parent
source = json.loads((HERE/'results/007-milp-search.json').read_text())
cases = [{'I': c['I'], 'J': c['J']} for c in source['result']['candidates']
         if c['domination_status'] == 'unknown']
cases += [{'I': [0, 1], 'J': [0, 2]},
          {'I': [0, 1, 2, 3, 5, 11], 'J': [0, 1, 3, 4, 5, 7, 13]}]
results = [find_certificate(c['I'], c['J']) for c in cases]
assert len(results) == 6
for c in results[:4]:
    verify(c)
assert all(c['status'] == 'unknown' for c in results[4:])

bad = []
c = deepcopy(results[0]); c['z'][0] = '99'; bad.append(c)
c = deepcopy(results[0]); c['variables'].reverse(); bad.append(c)
c = deepcopy(results[0]); c['terms'][0]['diagonals'] = [0, 0]; bad.append(c)
c = deepcopy(results[0])
for term in c['terms']: term['weight'] = str(-Fraction(term['weight']))
bad.append(c)
c = deepcopy(results[0])
for term in c['terms']: term['weight'] = '0'
bad.append(c)
for c in bad:
    try:
        verify(c)
    except AssertionError:
        pass
    else:
        raise AssertionError('Accepted corrupted certificate')

# An exact positive kernel proves incompleteness of this rule alone.
variables, rows = rows_for_support(cases[-1]['I'], cases[-1]['J'])
h = [2, 1, 1, 1, 1, 1, 3, 1, 1]
assert len(h) == len(variables) and all(t > 0 for t in h)
assert len(rows) == 3
assert all(sum(a*b for a, b in zip(row['row'], h)) == 0 for row in rows)
# Independently reconstruct all cancelled equalities at the rational point 2^-h.
point = {tuple(v): Fraction(1, 2**t) for v, t in zip(variables, h)}
I, J = cases[-1]['I'], cases[-1]['J']
def evaluated_terms(k):
    return sorted(point.get(('a', i), Fraction(1))*point.get(('b', k-i), Fraction(1))
                  for i in I if k-i in J)
for row in rows:
    k, ell = row['diagonals']
    assert sum(evaluated_terms(k)) == sum(evaluated_terms(ell))
assert sum(evaluated_terms(1)) == Fraction(3, 4)  # Not a solution of the full system.
print(json.dumps({'python': platform.python_version(), 'scipy': scipy.__version__,
                  'sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['monomial_cycles.py', 'verify_monomial_cycles.py',
                                          'check_monomial_cycles.py']},
                  'cases': results, 'adversarial_rejections': len(bad),
                  'positive_kernel': {'variables': variables, 'rows': rows, 'h': h,
                                      'full_diagonal_1': '3/4'}}, indent=2))
