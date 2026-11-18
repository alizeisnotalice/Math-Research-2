"""Recheck all nine-point residual certificates from a saved complete scan."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
COMPUTE = ROOT.parents[1]


def points(mask):
    return [i for i in range(mask.bit_length()) if (mask >> i) & 1]


def main():
    if not __debug__:
        raise RuntimeError('run without -O')
    production = ROOT / 'nine35-36-mac-return/summary.json'
    data = json.loads(production.read_text())
    residual = ROOT / 'nine-residual'
    report = json.loads((residual / 'summary.json').read_text())
    assert hashlib.sha256(production.read_bytes()).hexdigest() == report['input_sha256']
    expected = [(points(a), points(b)) for p in data['pairs'].values() for a, b in p['unknown_masks']]
    assert len(expected) == len({(tuple(a),tuple(b)) for a,b in expected}) == 166
    rows = report['rows']
    assert [(row['I'], row['J']) for row in rows] == expected
    module = importlib.util.spec_from_file_location('diagonal', COMPUTE / 'diagonal_domination.py')
    diagonal = importlib.util.module_from_spec(module)
    module.loader.exec_module(diagonal)
    direct = positive = 0
    for index, (row, (a, b)) in enumerate(zip(rows, expected)):
        assert row['i'] == index
        if row['status'] == 'diagonal_verified':
            assert diagonal.verify_certificate(a, b, row['certificate'])
            direct += 1
        else:
            assert row['status'] == 'positive_verified' and row['file'] == f'candidate-{index:03}.json'
            path = residual / row['file']
            raw = path.read_bytes()
            assert hashlib.sha256(raw).hexdigest() == row['sha256']
            case = json.loads(raw)['case']
            assert (case['I'], case['J']) == (a, b)
            verification = subprocess.check_output([sys.executable, str(COMPUTE / 'verify_positive_certificate.py'), str(path)], text=True)
            assert json.loads(verification)['status'] == 'verified'
            positive += 1
    assert (direct, positive) == (163, 3)
    print(json.dumps({'status':'verified', 'diagonal':direct, 'positive':positive, 'unresolved':0}))


if __name__ == '__main__':
    main()
