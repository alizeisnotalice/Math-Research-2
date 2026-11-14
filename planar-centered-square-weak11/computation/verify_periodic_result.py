#!/usr/bin/env python3
"""Exact saved-data check, with optional strict improvement over a seed enclosure."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    for field in ('input', 'primary', 'independent'):
        parser.add_argument('--' + field, type=Path, required=True)
    parser.add_argument('--baseline-primary', type=Path)
    args = parser.parse_args()
    data, p, q = read(args.input), read(args.primary), read(args.independent)
    px, py, lam = (F(data[k]) for k in ('Px', 'Py', 'lam'))
    mass = sum((F(row[2]) for row in data['pts']), F(0))
    radius = F(p['tail_radius'])
    tail = mass / (px * py) * (1 + px / (2 * radius)) * (1 + py / (2 * radius))
    assert min(px, py, lam, mass, radius) > 0
    assert tail == F(p['tail_bound']) == F(q['tail_bound_at_cap']['fraction']) < lam
    assert mass == F(p['mass']) == F(q['base_mass'])
    assert p['input_sha256'] == q['input_sha256'] == sha(args.input)
    directory = Path(__file__).resolve().parent
    assert p['code_sha256'] == sha(directory / 'periodic_halo.py')
    assert q['script_sha256'] == sha(directory / 'independent_periodic_halo.py')
    for side, bound in (('inner', 'lower'), ('outer', 'upper')):
        area, ratio = F(p[side + '_area']), F(p[side + '_ratio'])
        assert area == F(q['area_' + bound]['fraction'])
        assert ratio == F(q['weak_ratio_' + bound]['fraction']) == lam * area / mass
    lower, upper = F(p['inner_ratio']), F(p['outer_ratio'])
    assert lower <= upper
    output = {'status': 'passed', 'scope': 'saved exact arithmetic and file identity',
              'lower': str(lower), 'seed_upper': str(upper), 'tail_bound': str(tail)}
    if args.baseline_primary:
        base = read(args.baseline_primary)
        difference = lower - F(base['outer_ratio'])
        assert difference > 0
        output['gain_over_previous_seed_upper'] = str(difference)
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
