#!/usr/bin/env python3
"""Recompute the complete finite-template dual and reject corrupted variants."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

here = Path(__file__).resolve().parent
source = here/'results/010-fullquadratic-dual-720.json'
report = json.loads(source.read_text())
assert report['basis']['multiplier_total_degree_at_most'] == 1
assert report['basis']['positive_atom_total_degree_at_most'] == 2
verifier = here/'verify_template_dual.py'
run = subprocess.run([sys.executable, str(verifier), str(source)], capture_output=True, text=True, check=True)
verified = json.loads(run.stdout)
assert verified['annihilated_multiplier_columns'] == 308
assert verified['positive_atoms_checked'] == 377
assert verified['minimum_atom_value'] == [1, 1]
assert Path('/Volumes/data').is_mount()
temp = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round010/dual-attacks')
temp.mkdir(parents=True, exist_ok=True)
attacks = {}
r = deepcopy(report); r['functional'] = []; attacks['empty'] = r
r = deepcopy(report)
for item in r['functional']:
    item['coefficient'][0] *= -1
attacks['negative_functional'] = r
r = deepcopy(report)
for item in r['functional']:
    if item['monomial'] == []:
        item['coefficient'][0] += item['coefficient'][1]
        break
else:
    raise AssertionError('Expected nonzero constant coordinate')
attacks['broken_annihilation'] = r
r = deepcopy(report); r['unit_variables'].append('a2'); r['strict_variables'].remove('a2')
attacks['forged_unit_partition'] = r
r = deepcopy(report); r['functional'].append(deepcopy(r['functional'][0]))
attacks['duplicate_coordinate'] = r
for name, bad in attacks.items():
    path = temp/(name+'.json'); path.write_text(json.dumps(bad))
    run = subprocess.run([sys.executable, str(verifier), str(path)], capture_output=True, text=True)
    assert run.returncode != 0, name
squarefree = subprocess.run([sys.executable, str(verifier), str(here/'results/010-template-dual-720.json')], capture_output=True, text=True, check=True)
assert json.loads(squarefree.stdout)['positive_atoms_checked'] == 338
print(json.dumps({'verified_fullquadratic': verified, 'verified_squarefree': json.loads(squarefree.stdout), 'rejected_mutations': list(attacks),
                  'sha256': {f: hashlib.sha256((here/f).read_bytes()).hexdigest()
                             for f in ['template_dual_search.py', 'verify_template_dual.py',
                                       'check_template_dual.py']}}, indent=2))
