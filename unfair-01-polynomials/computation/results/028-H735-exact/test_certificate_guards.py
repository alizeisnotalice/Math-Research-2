"""Fast negative tests for the H7(3,5) independent certificate guards."""
import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from fractions import Fraction as F

from verify_frontend import verify, sha
from verify_h735_ball_certificate import (
    balls_from_boxes, coefficient, req, validate_certificate_metadata,
    verify_modular_galois,
)


def rejected(action):
    try:
        action()
    except (ValueError, KeyError, TypeError):
        return
    raise AssertionError('tampered input was accepted')


def main():
    root = Path(__file__).resolve().parent
    balls_path = root / 'H7-a3b5-root-balls.json'
    balls = json.loads(balls_path.read_text())
    front, _, boxes = verify(balls)
    metadata = {'balls_sha256': sha(balls_path),
                'canonical_sha256': front['canonical_sha256'], 'bits': 128}
    validate_certificate_metadata(metadata, balls_path, front['canonical_sha256'])
    bool_bits = copy.deepcopy(metadata); bool_bits['bits'] = True
    wrong_roots = copy.deepcopy(metadata); wrong_roots['balls_sha256'] = '0' * 64
    rejected(lambda: validate_certificate_metadata(bool_bits, balls_path, front['canonical_sha256']))
    rejected(lambda: validate_certificate_metadata(wrong_roots, balls_path, front['canonical_sha256']))

    # A positive constant coefficient cannot be silently marked as a negative
    # witness: the exact directed-ball guard has centre D and radius zero.
    D, factors = balls_from_boxes({0: [(F(1), F(1))]}, 11)
    centre, radius = coefficient(factors, [0], 0, D)
    assert (centre, radius) == (D, 0)
    rejected(lambda: req(centre + radius < 0, 'positive constant as negative witness'))

    galois_path = root / 'H7-mother-certificate.json'
    verify_modular_galois(galois_path)
    with TemporaryDirectory(dir=Path(__file__).parent) as temporary:
        tampered = json.loads(galois_path.read_text())
        tampered['certificates'][0]['coefficients_ascending'][3] = 1
        altered_path = Path(temporary) / 'wrong-H.json'
        altered_path.write_text(json.dumps(tampered))
        rejected(lambda: verify_modular_galois(altered_path))

    result = {
        'status': 'certificate_guard_negative_tests_passed',
        'tamper_rejections': [
            'bool_source_bits', 'wrong_root_box_sha', 'positive_constant_marked_negative',
            'wrong_H_in_galois_certificate'
        ],
    }
    Path(__file__).with_name('certificate-guard-tests.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
