"""Schema-only tests for the compact positional witness convention."""
import copy
import json
from pathlib import Path
from verify_frontend import verify
from verify_h735_ball_certificate import compact_records, standard_mask


def rejected(fn):
    try:
        fn()
    except (ValueError, KeyError, TypeError):
        return True
    raise AssertionError('invalid compact certificate accepted')


def main():
    root = Path(__file__).resolve().parent
    balls = json.loads((root / 'H7-a3b5-root-balls.json').read_text())
    _, layouts, _ = verify(balls)
    standard = standard_mask(balls)
    certificate = {'format': 'H7-a3b5-compact-ball-witness-v1', 'status': 'complete',
                   'witnesses': [[-1, -1] if mask == standard else [0, 0]
                                 for _, mask in layouts]}
    records = compact_records(certificate, layouts, standard)
    assert len(records) == 119062
    bad_bool = copy.deepcopy(certificate); bad_bool['witnesses'][0][0] = True
    bad_standard = copy.deepcopy(certificate)
    position = next(i for i, (_, mask) in enumerate(layouts) if mask == standard)
    bad_standard['witnesses'][position] = [0, 0]
    bad_length = copy.deepcopy(certificate); bad_length['witnesses'].pop()
    assert rejected(lambda: compact_records(bad_bool, layouts, standard))
    assert rejected(lambda: compact_records(bad_standard, layouts, standard))
    assert rejected(lambda: compact_records(bad_length, layouts, standard))
    result = {'status': 'compact_positional_schema_passed', 'witness_count': len(records),
              'tamper_rejections': ['bool', 'wrong_standard', 'wrong_length']}
    Path(__file__).with_name('compact-schema-tests.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
