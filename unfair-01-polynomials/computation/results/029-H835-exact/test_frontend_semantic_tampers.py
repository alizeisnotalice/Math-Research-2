import os
"""Direct semantic tamper regressions for H8(3,5), independent of hashes."""
import copy
import json
from pathlib import Path
import struct
from tempfile import TemporaryDirectory

from verify_frontend import verify


ROOT = Path(__file__).resolve().parent
BALLS = ROOT / 'H8-a3b5-root-balls.json'
LEDGER = Path(os.environ['H835_INPUT_DIR']) / 'H8-a3b5-canonical-masks.u32le'


def rejected(action, label):
    try:
        action()
    except (ValueError, KeyError, TypeError):
        return label
    raise AssertionError(f'accepted semantic tamper: {label}')


def main():
    original = json.loads(BALLS.read_text())
    duplicate = copy.deepcopy(original)
    # Groups 0 and 1 have the same mother and exponent.  Copy the Rouché
    # witness too, so rejection must use target-root disk separation.
    for key in ('disk', 'lhs', 'rhs'):
        duplicate['groups'][1][key] = copy.deepcopy(duplicate['groups'][0][key])
    wrong_block = copy.deepcopy(original)
    wrong_block['groups'][0]['block'] = 1
    wrong_kind = copy.deepcopy(original)
    wrong_kind['groups'][0]['kind'] = 'linear'
    rejected(lambda: verify(duplicate), 'duplicate_lift_disk_with_synced_rouche')
    rejected(lambda: verify(wrong_block), 'wrong_mother_child_binding')
    rejected(lambda: verify(wrong_kind), 'quadratic_group_relabelled_linear')

    data = bytearray(LEDGER.read_bytes())
    with TemporaryDirectory(dir=Path(os.environ['H835_INPUT_DIR'])) as temporary:
        directory = Path(temporary)
        duplicate_words = bytearray(data)
        duplicate_words[-4:] = duplicate_words[:4]
        duplicate_path = directory / 'duplicate.u32le'
        duplicate_path.write_bytes(duplicate_words)
        rejected(lambda: verify(original, duplicate_path), 'canonical_duplicate_and_omission')
        replacement_words = bytearray(data)
        replacement_words[-4:] = struct.pack('<I', 0)
        replacement_path = directory / 'replacement.u32le'
        replacement_path.write_bytes(replacement_words)
        rejected(lambda: verify(original, replacement_path), 'canonical_same_count_replacement')

    result = {'status': 'H835_frontend_semantic_tampers_passed',
              'rejections': ['duplicate_lift_disk_with_synced_rouche', 'wrong_mother_child_binding',
                             'quadratic_group_relabelled_linear', 'canonical_duplicate_and_omission',
                             'canonical_same_count_replacement'],
              'notes': 'Direct semantic validation; no SHA mismatch is used.'}
    Path(__file__).with_name('H835-frontend-semantic-tampers.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
