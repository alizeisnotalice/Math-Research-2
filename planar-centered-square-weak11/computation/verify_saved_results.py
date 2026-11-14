#!/usr/bin/env python3
"""Check saved R001 identities and source binding; this does not rerun geometry."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = json.loads((BASE / 'results/round001/primary.json').read_text())
    q = json.loads((BASE / 'results/round001/independent.json').read_text())
    source = BASE / 'inputs/legacy_r004_transfer_34_27.json'
    data = json.loads(source.read_text())
    mass = sum((F(row[2]) for row in data['pts']), F(0))
    px, py, level = (F(data[k]) for k in ('Px', 'Py', 'lam'))
    radius = F(p['tail_radius'])
    tail = mass / (px * py) * (1 + px / (2 * radius)) * (1 + py / (2 * radius))
    assert tail < level
    assert p['input_sha256'] == q['input_sha256'] == sha(source)
    assert p['code_sha256'] == sha(BASE / 'periodic_halo.py')
    assert q['script_sha256'] == sha(BASE / 'independent_periodic_halo.py')
    assert mass == F(p['mass']) == F(q['base_mass'])
    assert tail == F(p['tail_bound']) == F(q['tail_bound_at_cap']['fraction'])
    for side, bound in (('inner', 'lower'), ('outer', 'upper')):
        area = F(p[side + '_area'])
        ratio = F(p[side + '_ratio'])
        assert area == F(q['area_' + bound]['fraction'])
        assert ratio == F(q['weak_ratio_' + bound]['fraction']) == area * level / mass
    lower, upper = F(p['inner_ratio']), F(p['outer_ratio'])
    inherited = json.loads((BASE.parent / 'references/inherited_interval.json').read_text())
    assert lower == F(inherited['lower']) < upper
    assert p['counts']['boxes'] == q['counts']['all_span_admissible_boxes'] == 2329036
    assert all(q['primary_comparison']['checks'].values())
    assert all(q['primary_rectangle_comparison']['checks'].values())
    print(json.dumps({'status': 'passed', 'scope': 'saved exact arithmetic and file identity only',
                      'seed_ratio_gap': str(upper - lower),
                      'seed_upper_is_general_C2_upper': False}, indent=2))


if __name__ == '__main__':
    main()
