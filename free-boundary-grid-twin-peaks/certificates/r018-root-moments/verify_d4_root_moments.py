#!/usr/bin/env python3
"""Audit the D4-orbit moment accumulator against exact known certificates."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

KNOWN = {
    4: (779299174080, 5272309842240),
    5: (7063261481362958280000, 67631454448823894520000),
    6: (304351937651819767235871772704107520,
        3934566741149040698531461804200320000),
    7: (137077928808481869759712347602605837601908319657172480,
        2313395312707855007943501445611841392991226454185634816),
}

out = {'d4_runs': []}
for n, (a2_expected, d_expected) in KNOWN.items():
    x = json.loads(Path(f'n{n}-d4-root-moments.json').read_text())
    A, F, G, D = map(int, (x['A2'], x['firstsum'], x['secondsum'],
                            x['squared_euclidean_total']))
    assert A == a2_expected and D == d_expected
    assert F - 2*G == 4*D
    if n <= 6:
        y = json.loads(Path(f'n{n}-root-moments.json').read_text())
        assert (A, F, G) == tuple(map(int, (y['A2'], y['firstsum'], y['secondsum'])))
    if n == 7:
        assert x['threads'] == 8
    V, C = Fraction(F, 8*A), Fraction(G, 4*A)
    Q = Fraction(2, n*n)*(V-C)
    out['d4_runs'].append({'n': n, 'A2': str(A), 'D': str(D),
                           'F': str(F), 'G': str(G),
                           'V': str(V), 'C': str(C), 'Q': str(Q)})
for p in (Path('n7_root_moments_d4.cpp'),
          Path('../../computation/d4_peak_dp_generic.cpp')):
    out.setdefault('sha256', {})[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
Path('d4-root-moments-verified.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
