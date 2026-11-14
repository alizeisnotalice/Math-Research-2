#!/usr/bin/env python3
"""Build the actually used finite-difference proposal table; not derivatives."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    keyed = {}
    for row in json.loads(args.results.read_text()):
        t = row['transform']
        key = (t['kind'], t.get('i', -1), t.get('axis', -1), t.get('anchor', -1))
        keyed.setdefault(key, {})[(t['eps'], t['sign'])] = F(row['inner_ratio'])
    keys = sorted(keyed)
    eps = ['1/1000', '1/10000', '1/100000']
    gradients = [[float((keyed[k][e, 1] - keyed[k][e, -1]) / (2 * F(e)))
                  for k in keys] for e in eps]
    data = {'schema': 'finite-difference-proposals-only', 'keys': [list(k) for k in keys],
            'eps': eps, 'gradients': gradients, 'budgets': [1, 2, 4, 8, 16, 32, 64, 154],
            'groups': ['all', 'mass', 'coordinate', 'geometry'],
            'trust_scales': ['1/1000000', '3/1000000', '1/100000', '3/100000',
                             '1/10000', '3/10000', '1/1000']}
    data['kinds'] = [row[0] for row in data['keys']]
    args.output.write_text(json.dumps(data, indent=2) + '\n')


if __name__ == '__main__':
    main()
