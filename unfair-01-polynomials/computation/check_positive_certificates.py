"""Black-box negative tests for literature_audit's exact certificate verifier."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import hashlib
import platform
import scipy
import sympy
from positive_certificate_search import discover, encode_monomial

base = Path(__file__).resolve().parent
verifier = base / 'verify_positive_certificate.py'
source = base / 'results/010-positive-719b.json'
assert Path('/Volumes/data').is_mount()
out = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round010/attacks')
out.mkdir(parents=True, exist_ok=True)
valid = json.loads(source.read_text())

def mutate_unit(r):
    r['unit_variables'].append('a2')
    r['strict_variables'].remove('a2')

def omit_strict(r):
    r['strict_variables'].remove('a2')

def nonzero_minimum(r):
    r['case']['I'] = [x+1 for x in r['case']['I']]

def nonstrict_multiplier(r):
    r['certificate']['multipliers'][0]['monomial'].append('a0')

def squarefree_atom_violation(r):
    atom = r['certificate']['positive_atoms'][0]
    atom['positive'] = ['a2', 'a2']

def atom_overlap(r):
    atom = r['certificate']['positive_atoms'][0]
    atom['positive'] = ['a2']
    atom['one_minus'] = ['a2']

def negative_atom_weight(r):
    r['certificate']['positive_atoms'][0]['coefficient'] = [-1, 1]

def false_identity(r):
    coefficient = r['certificate']['multipliers'][0]['coefficient']
    coefficient[0] += coefficient[1]

attacks = {
    'forged_unit': mutate_unit,
    'omitted_strict_variable': omit_strict,
    'support_not_normalized': nonzero_minimum,
    'nonstrict_multiplier': nonstrict_multiplier,
    'repeated_atom_variable': squarefree_atom_violation,
    'overlapping_atom_factors': atom_overlap,
    'negative_atom_weight': negative_atom_weight,
    'false_polynomial_identity': false_identity,
}

valid_path = out / 'valid.json'
valid_path.write_text(json.dumps(valid))
assert subprocess.run([sys.executable, str(verifier), str(valid_path)], capture_output=True, text=True).returncode == 0
result = {}
for name, attack in attacks.items():
    record = copy.deepcopy(valid)
    attack(record)
    path = out / f'{name}.json'
    path.write_text(json.dumps(record))
    run = subprocess.run([sys.executable, str(verifier), str(path)], capture_output=True, text=True)
    result[name] = {'returncode': run.returncode, 'rejected': run.returncode != 0,
                    'message': (run.stderr or run.stdout).strip()[:180]}
assert all(item['rejected'] for item in result.values()), result
positive_results = {}
for name in ['010-positive-718.json', '010-positive-719b.json']:
    record = json.loads((base/'results'/name).read_text())
    assert all(len(t['monomial']) <= 1 for t in record['certificate']['multipliers'])
    assert all(len(t['positive'])+len(t['one_minus']) <= 2
               for t in record['certificate']['positive_atoms'])
    run = subprocess.run([sys.executable, str(verifier), str(base/'results'/name)],
                         capture_output=True, text=True, check=True)
    positive_results[name] = json.loads(run.stdout)
assert discover([0, 1], [0, 1], 1, 2, 4096)['status'] == 'unit_collision_infeasible'
assert discover([0, 1], [0, 2], 1, 2, 4096)['status'] == 'trivially_uniform_no_strict_variables'
assert discover([0], [0], 1, 2, 4096)['status'] == 'trivially_uniform_no_strict_variables'
assert encode_monomial(sympy.Symbol('a2')**2) == ['a2', 'a2']
print(json.dumps({'valid': 'accepted', 'attacks': result, 'verified': positive_results,
                  'unit_boundary_tests': 3, 'square_multiplier_encoding': True,
                  'python': platform.python_version(), 'sympy': sympy.__version__,
                  'scipy': scipy.__version__,
                  'sha256': {f: hashlib.sha256((base/f).read_bytes()).hexdigest()
                             for f in ['positive_certificate_search.py', 'verify_positive_certificate.py',
                                       'check_positive_certificates.py']}}, indent=2))
