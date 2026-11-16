#!/usr/bin/env python3
"""Portable independent n=3..5 reference check for this archived package."""
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
REF = HERE.parent / 'fixed-state-completion-root-r023' / 'top2-reference'

def load(path):
    return json.loads(path.read_text())

def compare(n, source_output):
    d4 = load(HERE / f'n{n}.json')
    archived = load(REF / f'top2-reference-n{n}.json')
    for key, other in [('Btop2', 'Btop2'),
                       ('tau2_squared_euclidean_total', 'D'),
                       ('tau2_firstsum', 'F'),
                       ('tau2_secondsum', 'G')]:
        assert int(d4[key]) == int(archived[other]) == int(source_output[other]), (n, key)
    h = [int(v) for v in source_output['hist']]
    assert h == [int(d4['tau2_distance_counts'].get(str(i), 0)) for i in range(2 * n - 1)]
    assert int(source_output['F']) - 2 * int(source_output['G']) == 4 * int(source_output['D'])
    return hashlib.sha256((REF / f'top2-reference-n{n}.json').read_bytes()).hexdigest()

def check_summary():
    summary = load(HERE / 'summary-n3-n7.json')
    assert [x['n'] for x in summary['results']] == [3, 4, 5, 6, 7]
    for x in summary['results']:
        n = x['n']; raw = load(HERE / f'n{n}.json')
        q = Fraction(int(raw['tau2_squared_euclidean_total']), n * n * int(raw['Btop2']))
        got = x['tau2_Q']
        assert q == Fraction(int(got['numerator']), int(got['denominator'])), n

def main():
    with tempfile.TemporaryDirectory(dir='/Volumes/data/CodexTemp/T2') as td:
        exe = pathlib.Path(td) / 'top2_subset_reference'
        subprocess.run(['c++', '-O3', '-std=c++20', str(REF / 'top2_subset_reference.cpp'), '-o', str(exe)], check=True)
        hashes = []
        for n in (3, 4, 5):
            actual = json.loads(subprocess.check_output([str(exe), str(n)], text=True))
            hashes.append(compare(n, actual))
    check_summary()
    print(json.dumps({'passed': True, 'independent_n': [3, 4, 5],
                      'corrected_Q_summary': 'summary-n3-n7.json',
                      'reference_json_sha256': hashes}, indent=2))

if __name__ == '__main__':
    main()
