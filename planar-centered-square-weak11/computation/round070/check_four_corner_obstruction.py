#!/usr/bin/env python3
"""Exact finite-domain check; no halo or envelope quadrature."""
from fractions import Fraction as F
import hashlib
import json
import platform
from pathlib import Path
import sys
import time


def main():
    started = time.perf_counter()
    d = F(99, 100)
    points = [(F(0), F(0), F(1, 4)), (d, F(0), F(1, 4)),
              (F(0), d, F(1, 4)), (d, d, F(1, 4)),
              (F(33, 32), F(99, 200), F(1, 1000))]
    eligible = []
    rows = []
    for mask in range(1, 32):
        ps = [p for i, p in enumerate(points) if mask >> i & 1]
        mass = sum(p[2] for p in ps)
        span = max(max(p[j] for p in ps)-min(p[j] for p in ps) for j in (0, 1))
        ok = span * span < mass
        rows.append({'mask': mask, 'mass': str(mass), 'span_squared': str(span*span), 'eligible': ok})
        if ok:
            eligible.append((mass, mask))
    best = max(m for m, _ in eligible)
    maximizers = [mask for m, mask in eligible if m == best]
    assert best == 1 and maximizers == [15]
    reached = 1
    while True:
        nxt = reached
        for _, mask in eligible:
            if mask & reached:
                nxt |= mask
        if nxt == reached:
            break
        reached = nxt
    assert reached == 31
    checks = {
        'radical_square_margin': 7*F(199,50)**2-F(21,2)**2,
        'whole_nonterminal_margin': F(33,32)**2-F(35,33)*F(1001,1000),
        'whole_ineligible_margin': F(33,32)**2-F(1001,1000),
        'connector_pair_eligibility_margin': F(251,1000)-F(99,200)**2,
        'fat_margin': d-F(3,4),
        'support_tree_failure_margin': 3*d-F(5,2),
    }
    assert all(v > 0 for v in checks.values())
    assert checks['radical_square_margin'] == F(791,1250)
    result = {
        'status': 'PASS', 'python': platform.python_version(),
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input': [[str(v) for v in p] for p in points],
        'subset_count': 31, 'eligible_count': len(eligible),
        'unique_maximum_mask': 15, 'maximum_mass': '1',
        'connected': True, 'strict_rational_margins': {k: str(v) for k,v in checks.items()},
        'subsets': rows, 'elapsed_seconds': time.perf_counter()-started,
        'scope': 'Exact domain and scalar checks only; no halo areas or full envelope areas computed.'
    }
    out = json.dumps(result, indent=2)+'\n'
    if len(sys.argv)>1:
        Path(sys.argv[1]).write_text(out)
    else:
        print(out, end='')


if __name__ == '__main__':
    main()
