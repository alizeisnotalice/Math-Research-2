#!/usr/bin/env python3
"""Exact audit of centered vertical-cut covariance accumulators."""
import json
from pathlib import Path

EXPECTED_A2 = {3: 143808, 4: 779299174080, 5: 7063261481362958280000}
pair = {z['n']: z for z in json.loads(Path('peak-association-exact.json').read_text())['results']}
out = {'runs': []}
for n, expected in EXPECTED_A2.items():
    z = json.loads(Path(f'n{n}-cut-ie.json').read_text())
    A = int(z['A2']); assert A == expected
    m = list(map(int, z['cut_mean_numerators']))
    T = [[int(x) for x in row] for row in z['cut_ordered_pair_numerators']]
    assert len(m) == n-1 and len(T) == n-1 and all(len(r) == n-1 for r in T)
    cov = [[2*A*T[j][k]-m[j]*m[k] for k in range(n-1)] for j in range(n-1)]
    assert all(cov[j][k] == cov[k][j] for j in range(n-1) for k in range(n-1))
    assert all(x < 0 for row in cov for x in row)
    if n in pair:
        old = pair[n]['centered_cut_covariances']
        assert cov == [[x['numerator_over_4A2'] for x in row] for row in old]
    out['runs'].append({'n': n, 'A2': str(A), 'cut_mean_numerators': list(map(str,m)),
                        'centered_covariance_numerators_over_4A2': [[str(x) for x in row] for row in cov],
                        'all_entries_strictly_negative': True,
                        'n3_n4_pair_dp_crosscheck': n in pair})
Path('cut-covariances-verified.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
