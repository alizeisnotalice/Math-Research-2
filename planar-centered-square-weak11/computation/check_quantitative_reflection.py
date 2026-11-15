#!/usr/bin/env python3
"""Exact rational checks for quantitative-reflection-gap.md; no numerical halo claim."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import platform
import sys

assert os.path.ismount('/Volumes/data'), 'data volume is not mounted'
HERE = Path(__file__).resolve()
TOPIC = HERE.parents[1]
ROOT = TOPIC.parent
OUTPUT = TOPIC / 'computation/results/quantitative-reflection/exact-checks.json'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rat(x):
    return f'{x.numerator}/{x.denominator}'

checks = {
    'k_squared_half_exceeds_one_twelfth': F(17,6)**2 > 8,
    'k_squared_quarter_exceeds_one_twentyfifth': F(71,25)**2 > 8,
}
d = F(1,12)
rows = []
for n in range(1,9):
    closed = F(2160, 25920**(2**(n-1)))
    assert d == closed
    t = d / 45
    tau = d*d/180
    next_d = d*d/2160
    row_checks = {
        'positive_and_at_most_nine': 0 < d <= 9,
        't_positive_and_at_most_one_fifth': 0 < t <= F(1,5),
        'large_mass_branch': tau/F(12) == next_d,
        'small_mass_first_overhead': 9*t/(1-t) <= d/4,
        'small_mass_second_overhead': tau/t == d/4,
        'recursive_gap_at_most_half': next_d <= d/2,
        'closed_form_equal': closed == d,
    }
    assert all(row_checks.values())
    rows.append({'N': n, 'gap': rat(d), 't': rat(t), 'tau': rat(tau), 'checks': row_checks})
    d = next_d
assert all(checks.values())
source_paths = [
    TOPIC/'proofs/reflection-nonattainment.md',
    Path('/Volumes/data/CodexProjects/T/planar-centered-square-weak11/scratch/round005-public-stage/round005/proof.md'),
]
result = {
    'status': 'PASS',
    'python': sys.version,
    'platform': platform.platform(),
    'arithmetic': 'fractions.Fraction exact rational arithmetic; positive square comparisons for sqrt(2)',
    'scope': 'Scalar constants and N=1..8 recurrence identities only; general geometry and induction are analytic proofs.',
    'seed': None, 'floating_tolerance': None,
    'command': f'python3 {HERE}',
    'script_sha256': sha(HERE),
    'proof_sha256': sha(TOPIC/'proofs/quantitative-reflection-gap.md'),
    'sources': [{'path': str(p), 'sha256': sha(p)} for p in source_paths],
    'scalar_checks': checks, 'recurrence_checks': rows,
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':result['status'],'output':str(OUTPUT),'N':list(range(1,9)),'gap_N2':rows[1]['gap'],'script_sha256':result['script_sha256']},ensure_ascii=False))
