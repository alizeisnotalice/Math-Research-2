#!/usr/bin/env python3
"""Check R002 saved exact improvements and input changes, not recompute geometry."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    original = read(BASE / 'inputs/legacy_r004_transfer_34_27.json')
    filename = BASE / 'inputs/r002_final_0475.json'
    final = read(filename)
    assert all(F(original[k]) == F(final[k]) for k in ('Px', 'Py', 'lam'))
    assert original['pts'][:50] == final['pts'][:50]
    assert all(F(p[2]) == F(q[2]) for p, q in zip(original['pts'], final['pts']))
    assert len(original['pts']) == len(final['pts']) == 51
    assert all(0 <= F(p[axis]) < F(final[key])
               for p in final['pts'] for axis, key in ((0, 'Px'), (1, 'Py')))
    result = BASE / 'results/round002'
    p = read(result / 'primary_final.json')
    q = read(result / 'independent_final.json')
    assert p['input_sha256'] == q['input_sha256'] == sha(filename)
    assert p['code_sha256'] == sha(BASE / 'periodic_halo.py')
    assert q['script_sha256'] == sha(BASE / 'independent_periodic_halo.py')
    assert F(p['tail_bound']) < F(final['lam'])
    for side, bound in (('inner', 'lower'), ('outer', 'upper')):
        assert F(p[side + '_area']) == F(q['area_' + bound]['fraction'])
        assert F(p[side + '_ratio']) == F(q['weak_ratio_' + bound]['fraction'])
    lower = F(p['inner_ratio'])
    old = F(read(BASE.parent / 'references/inherited_interval.json')['lower'])
    previous_upper = F(read(result / 'independent_best.json')['weak_ratio_upper']['fraction'])
    assert lower > previous_upper > old
    ledger = read(result / 'search_ledger.json')
    assert sum(x['count'] for x in ledger) == 6305
    assert all(x['errors'] == 0 and x['workers'] == 12 for x in ledger)
    assert ledger[0]['scan_code_sha256'] == sha(BASE / 'scan_perturbations.py')
    assert all(x['scan_code_sha256'] == sha(BASE / 'scan_position_grid.py') for x in ledger[1:])
    print(json.dumps({'status': 'passed', 'lower': str(lower),
                      'gain_over_legacy_L4': str(lower - old),
                      'gain_over_intermediate_upper': str(lower - previous_upper),
                      'global_optimality': 'not proved'}, indent=2))


if __name__ == '__main__':
    main()
