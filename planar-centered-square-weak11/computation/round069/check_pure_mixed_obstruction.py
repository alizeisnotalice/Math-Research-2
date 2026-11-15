#!/usr/bin/env python3
"""Exact finite-domain and entire-patch checks; no complete halo area run."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import platform
import sys
import time


def main():
    assert __debug__ and os.path.ismount('/Volumes/data')
    started = time.perf_counter()
    points = [(F(0), F(0)), (F(99, 100), F(99, 100)),
              (F(0), F(99, 100)), (F(99, 100), F(0)),
              (F(-3, 10), F(0))]
    masses = [F(2499, 5000), F(1, 2), F(1, 10000), F(1, 10000), F(1, 2)]
    S = 15
    rows = []
    for mask in range(1, 32):
        ids = [i for i in range(5) if mask >> i & 1]
        mass = sum((masses[i] for i in ids), F(0))
        spans = [max(points[i][j] for i in ids)-min(points[i][j] for i in ids) for j in (0, 1)]
        rows.append({'mask': mask, 'mass': str(mass), 'span_squared': str(max(spans)**2), 'eligible': mass > max(spans)**2})
    top = max(F(r['mass']) for r in rows if r['eligible'])
    maxima = [r['mask'] for r in rows if r['eligible'] and F(r['mass']) == top]
    assert top == 1 and maxima == [S]
    assert rows[S-1]['eligible'] and rows[16]['eligible']  # fullS and {0,4} connect all sources
    full_mass = sum(masses)
    nonterminal_margin = F(129, 100)**2-F(35, 33)*full_mass
    assert nonterminal_margin > 0 and F(99, 100) > F(3, 4)
    side_cost = F(99, 100)
    diagonal_cost = 2*side_cost-side_cost**2
    assert diagonal_cost > side_cost and 3*side_cost > F(5, 2)

    radius = F(4999, 10000)
    xlo, xhi = -radius, radius-F(3, 10)
    ycut = F(177, 500)
    assert 0 < ycut < radius
    captured = masses[0]+masses[4]
    # Every open patch point captures {origin,external} at this same radius.
    full_witness_margin = captured-4*radius**2
    assert full_witness_margin > 0
    for i in (0, 4):
        assert max(abs(points[i][0]-xlo), abs(points[i][0]-xhi)) <= radius
    # No residual witness: vertical distance alone costs > total residual mass.
    residual_exclusion = 4*ycut**2-masses[4]
    origin_exclusion = 4*ycut**2-masses[0]
    assert residual_exclusion > 0 and origin_exclusion > 0
    # No S-witness containing a right-column atom, even using all S mass.
    right_gap = F(99, 100)-xhi
    right_exclusion = 4*right_gap**2-1
    assert right_exclusion > 0
    # After excluding the right column, all remaining S mass is left-column mass.
    top_gap = F(99, 100)-radius
    left_mass = masses[0]+masses[2]
    top_exclusion = 4*top_gap**2-left_mass
    assert top_exclusion > 0
    area = 2*(xhi-xlo)*(radius-ycut)
    assert area == F(5105041, 25000000) and area > F(1, 5)
    spec = {'id': 'r069-five-source-pure-mixed-over-one-fifth',
            'points': [list(map(str, p)) for p in points], 'initial_masses': list(map(str, masses)),
            'S_mask': S, 'threshold': '1', 'geometry': 'centered axis-parallel common-radius squares',
            'patches': [[str(xlo), str(xhi), str(ycut), str(radius)],
                        [str(xlo), str(xhi), str(-radius), str(-ycut)]],
            'common_half_side': str(radius), 'patch_endpoints': 'open'}
    root = Path(__file__).resolve().parent
    ip = root/'r069-pure-mixed-obstruction-input.json'
    ip.write_text(json.dumps(spec, indent=2)+'\n')
    result = {'schema': 'r069-pure-mixed-patch-arithmetic-v1', 'status': 'PASS',
              'python': sys.version, 'platform': platform.platform(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'input_sha256': hashlib.sha256(ip.read_bytes()).hexdigest(),
              'complete_subset_count': len(rows), 'subsets': rows, 'maximum_mass': str(top),
              'all_maximum_masks': maxima, 'full_input_mass': str(full_mass),
              'nonterminal_margin': str(nonterminal_margin), 'MST_cost': str(3*side_cost),
              'minimum_patch_witness_slack': str(full_witness_margin),
              'residual_exclusion_slack': str(residual_exclusion),
              'origin_exclusion_slack': str(origin_exclusion),
              'right_column_exclusion_slack': str(right_exclusion),
              'top_left_exclusion_slack': str(top_exclusion),
              'patch_area': str(area), 'excess_over_one_fifth': str(area-F(1, 5)),
              'wall_seconds': time.perf_counter()-started,
              'scope': 'Whole open rational patch lies in H(A) minus (H(S) union H(A\\S)); it strictly refutes |pure-mixed|<=L/5. No total halo area, net interaction credit bound or general3.5/C2 claim.'}
    (root/'r069-pure-mixed-obstruction-arithmetic.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'patch_area': str(area), 'excess': str(area-F(1, 5)), 'seconds': result['wall_seconds']}))


if __name__ == '__main__':
    main()
