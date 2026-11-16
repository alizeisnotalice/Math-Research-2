#!/usr/bin/env python3
"""Exact rational audit for independent-set root moment accumulators."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

EXPECTED = {
    4: (779299174080, 5272309842240),
    5: (7063261481362958280000, 67631454448823894520000),
    6: (304351937651819767235871772704107520,
        3934566741149040698531461804200320000),
}

out = {'runs': []}
for n, (expected_a2, expected_d) in EXPECTED.items():
    data = json.loads(Path(f'n{n}-root-moments.json').read_text())
    a2 = int(data['A2']); first = int(data['firstsum']); second = int(data['secondsum'])
    marg = list(map(int, data['root_marginals']))
    assert a2 == expected_a2
    assert sum(marg) == 2*a2
    assert all(x > 0 for x in marg)
    d = (first - 2*second) // 4
    assert first - 2*second == 4*d
    assert d == expected_d
    V = Fraction(first, 8*a2)
    C = Fraction(second, 4*a2)
    q_from_covariance = Fraction(2, n*n) * (V-C)
    q_from_distance = Fraction(d, n*n*a2)
    assert q_from_covariance == q_from_distance
    out['runs'].append({
        'n': n, 'independent_sets': data['independent_sets'], 'A2': str(a2),
        'firstsum': str(first), 'secondsum': str(second), 'D': str(d),
        'V_E_norm_X_minus_c_sq': str(V), 'C_E_dot_centered_roots': str(C),
        'Q': str(q_from_covariance),
        'root_marginal_sum': str(sum(marg)),
        'center_marginal': str(marg[(n*n)//2]) if n % 2 else None,
        'corner_marginal': str(marg[0]),
    })
for name in ('root_moments_ie.cpp',
             '../../computation/independent_peak_dp.cpp'):
    p = Path(name)
    out.setdefault('sha256', {})[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
Path('root-moments-verified.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
