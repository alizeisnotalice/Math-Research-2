#!/usr/bin/env python3
"""Integer checks for the n=3,4,5 root-depth ledgers."""
import json
from fractions import Fraction
from pathlib import Path

A2 = {3: 143808, 4: 779299174080, 5: 7063261481362958280000}
Q = {
    3: Fraction(1100, 2247),
    4: Fraction(784569917, 1855474224),
    5: Fraction(563595453740199121, 1471512808617282975),
}
rows=[]
for n,a2 in A2.items():
    d=json.loads(Path(f'root-depth-n{n}.json').read_text())
    h={int(k):int(v) for k,v in d['all_roots_histogram'].items()}
    assert sum(h.values()) == int(d['all_roots_total']) == 2*a2
    z=sum(k*k*v for k,v in h.items())
    assert z == int(d['all_roots_depth_squared_sum'])
    lower=Fraction(z,2*n*n*a2)
    assert 0 < lower <= Q[n]
    rows.append({'n':n,'depth_squared_lower_Q':str(lower),'fraction_of_actual_Q':str(lower/Q[n])})
out={'status':'passed','rows':rows}
Path('verify-root-depth.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
