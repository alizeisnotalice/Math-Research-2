"""Export exact finite-difference data for heuristic MATLAB LP proposals."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    raw = Path(args.input).read_bytes()
    data = json.loads(raw)
    columns = [[Fraction(x) for x in c['unit_mass_increment']]
               for c in data['columns']]
    assert all(sum(c) == 0 for c in columns)
    n = len(data['target_masses'])
    assert all(len(c) == n for c in columns)
    labels = ['1/1000', '1/10000']
    output = {
        'schema': 'mass-lp-numeric-proposals-v1',
        'scope': 'Heuristic finite differences, not certified derivatives or bounds',
        'exact_input_sha256': hashlib.sha256(raw).hexdigest(),
        'W': [[float(c[i]) for c in columns] for i in range(n)],
        'gradient_scales': labels,
        'G': [[float(Fraction(c['estimates'][s]['gradient_fraction_of_interval_midpoints']))
               for c in data['columns']] for s in labels],
        'mass': [float(Fraction(x)) for x in data['target_masses']],
        'lower': [float(Fraction(x)) for x in data['lower_mass_bound']],
        'upper': [float(Fraction(x)) for x in data['upper_mass_bound']],
        'budgets': data['budgets'],
        'trust_exact': data['trust_scales'],
        'trust': [float(Fraction(x)) for x in data['trust_scales']],
    }
    Path(args.output).write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
