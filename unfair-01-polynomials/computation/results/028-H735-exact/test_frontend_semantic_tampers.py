"""Semantic tamper regressions for the independent H7(3,5) front end.

These call the front-end verifier directly, so no SHA binding is involved.
"""
import copy
import gzip
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from verify_frontend import verify


ROOT = Path(__file__).resolve().parent
BALLS = ROOT / 'H7-a3b5-root-balls.json'
LEDGER = ROOT / 'H7-a3b5-canonical-masks.jsonl.gz'


def rejected(action, label):
    try:
        action()
    except (ValueError, KeyError, TypeError):
        return label
    raise AssertionError(f'semantic tamper accepted: {label}')


def invoke_with_ledger(rows, directory, label):
    ledger = directory / f'{label}.jsonl'
    ledger.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows))
    run = subprocess.run(
        [sys.executable, str(Path(__file__).with_name('verify_frontend.py')),
         '--balls', str(BALLS), '--record', str(directory / f'{label}.json'),
         '--masks', str(directory / f'{label}-masks.json'), '--source-ledger', str(ledger)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if run.returncode == 0:
        raise AssertionError(f'canonical ledger tamper accepted: {label}')
    return label


def main():
    original = json.loads(BALLS.read_text())

    # Indices 0 and 1 are both a-lift quadratic groups over mother block 0.
    # Copy its disk and its matching lift Rouché bounds together, so the
    # per-group Rouché calculation still succeeds; only global root-disk
    # separation may reject this forged repeated target root.
    duplicate = copy.deepcopy(original)
    first, second = duplicate['groups'][0], duplicate['groups'][1]
    for key in ('disk', 'lhs', 'rhs'):
        second[key] = copy.deepcopy(first[key])

    wrong_parent = copy.deepcopy(original)
    wrong_parent['groups'][0]['block'] = 1

    wrong_linear_kind = copy.deepcopy(original)
    linear = next(group for group in wrong_linear_kind['groups'] if group['index'] == 24)
    linear['kind'] = 'quadratic'

    rejected(lambda: verify(duplicate), 'duplicate_disk_with_synchronised_rouche')
    rejected(lambda: verify(wrong_parent), 'wrong_mother_child_block_binding')
    rejected(lambda: verify(wrong_linear_kind), 'linear_group_kind_changed')

    with gzip.open(LEDGER, 'rt') as handle:
        rows = [json.loads(line) for line in handle]
    assert len(rows) == 119062
    duplicate_row = copy.deepcopy(rows)
    duplicate_row[-1] = copy.deepcopy(duplicate_row[0])
    replacement = copy.deepcopy(rows)
    replacement[-1] = {'E': 0, 'mask': '0' * 29, 'selected_group_indices': []}
    with TemporaryDirectory(dir=Path(__file__).parent) as temporary:
        directory = Path(temporary)
        invoke_with_ledger(duplicate_row, directory, 'duplicate_omission')
        invoke_with_ledger(replacement, directory, 'same_count_replacement')

    result = {
        'status': 'frontend_semantic_tamper_regressions_passed',
        'rejections': [
            'duplicate_disk_with_synchronised_rouche',
            'wrong_mother_child_block_binding',
            'linear_group_kind_changed',
            'canonical_duplicate_and_omission',
            'canonical_same_count_replacement',
        ],
        'notes': 'All tests call semantic validators directly; none rely on a SHA mismatch.',
    }
    Path(__file__).with_name('frontend-semantic-tamper-tests.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
