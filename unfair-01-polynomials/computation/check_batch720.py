#!/usr/bin/env python3
"""Validate all returned 7x20 supports and their exact exclusion evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from diagonal_domination import verify_certificate

here = Path(__file__).resolve().parent
batch = json.loads((here/'results/010-structural-search720.json').read_text())
proofs = {}
def canonical(I, J):
    return min((tuple(I), tuple(J)),
               (tuple(sorted(I[-1]-i for i in I)), tuple(sorted(J[-1]-j for j in J))))
for i in range(1,5):
    path = here/f'results/010-positive-720r{i}.json'
    record = json.loads(path.read_text())
    subprocess.run([sys.executable, str(here/'verify_positive_certificate.py'), str(path)],
                   capture_output=True, text=True, check=True)
    proofs[canonical(record['case']['I'], record['case']['J'])] = path.name
run = subprocess.run([sys.executable, str(here/'check_escape720_identity.py')],
                     capture_output=True, text=True, check=True)
identity = json.loads(run.stdout)
proofs[canonical(identity['I'], identity['J'])] = '010-escape720-identity.json'
counts = Counter(); checked = []; seen = set()
for record in batch['result']['records']:
    if 'I' not in record:
        continue
    I, J = record['I'], record['J']
    assert I == sorted(set(I)) and J == sorted(set(J))
    assert I[0] == J[0] == 0 and I[-1] == 7 and J[-1] == 20
    assert len(I)+len(J) >= 11 and (set(I)&set(J))-{0}
    key = (tuple(I),tuple(J)); assert key not in seen; seen.add(key)
    multiplicity = Counter(i+j for i in I for j in J)
    assert sorted(k for k, v in multiplicity.items() if v == 1) == [0,7,20,27]
    if record['verdict']['status'] == 'diagonal_domination':
        verify_certificate(I, J, record['verdict']['certificate'])
        evidence = 'stored direct domination certificate'
        counts['direct_domination'] += 1
    else:
        evidence = proofs[canonical(I,J)]
        counts['positive_polynomial_or_bezout_and_reversal'] += 1
    checked.append({'I':I,'J':J,'exact_exclusion':evidence})
assert counts == {'direct_domination':36,'positive_polynomial_or_bezout_and_reversal':10}
print(json.dumps({'counts':dict(counts), 'distinct_returned_supports':len(seen),
                  'complete_support_enumeration_claimed':False, 'checked':checked,
                  'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
