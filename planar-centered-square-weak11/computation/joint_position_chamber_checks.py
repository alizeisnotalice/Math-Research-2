#!/usr/bin/env python3
"""Exact small checks for pure-position chambers and envelope error order.
Standard-library Fraction only; no halo search, optimizer or real candidate.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

start = time.process_time()
def text(q):
    return str(F(q))
def area(rect):
    x0, x1, y0, y1 = rect
    return max(F(0), x1 - x0) * max(F(0), y1 - y0)
def halo_rect(s):
    # Unique relevant periodic label for Px=Py=10, M=4, lambda=1, R=2.
    x, y, radius = F(5) + s, F(5), F(1)
    return x - radius, x + radius, y - radius, y + radius

tail = F(4, 100) * (1 + F(10, 4)) ** 2
assert tail == F(49, 100) < 1
rows = []
for n in range(1, 11):
    h = F(1, 1 << n)
    for s in [-h, F(0), h]:
        assert area(halo_rect(s)) == 4
    # The continuous family consists of equal-height translated rectangles.
    # Their exact union spans the extreme endpoints and has no x gaps.
    left, right = halo_rect(-h), halo_rect(h)
    union_area = area((left[0], right[1], left[2], left[3]))
    assert union_area == 4 + 4 * h
    ratio_upper, true_ratio = union_area / 4, F(1)
    assert (ratio_upper - true_ratio) / h == 1
    rows.append({'h': text(h), 'every_point_area_identity': '4',
                 'family_union_area': text(union_area),
                 'true_max_ratio': '1', 'family_union_ratio': text(ratio_upper),
                 'gap': text(ratio_upper - true_ratio), 'gap_over_h': '1'})

H = [[F(0), F(1)], [F(1), F(0)]]
def mv(m, v):
    return [sum((a * b for a, b in zip(row, v)), F(0)) for row in m]
def value(v):
    return 1 + v[0] * v[1]
for v, eigenvalue in [([F(1), F(1)], F(1)),
                      ([F(1), F(-1)], F(-1))]:
    assert mv(H, v) == [eigenvalue * a for a in v]
positive, negative = [F(3, 5), F(4, 5)], [F(3, 5), F(-4, 5)]
assert sum(a * a for a in positive) == sum(a * a for a in negative) == 1
assert value(positive) == F(37, 25) and value(negative) == F(13, 25)
assert value([F(21, 29), F(20, 29)]) - 1 == F(420, 841) < F(1, 2)
# Root-free algebra proves the spectral value: 2xy <= x^2+y^2,
# equality on x=y; on the unit circle x^2=y^2=1/2, hence gain=1/2.
for x, y in [(F(3, 5), F(4, 5)), (F(21, 29), F(20, 29)),
             (F(-3, 5), F(-4, 5))]:
    assert x * x + y * y - 2 * x * y == (x - y) ** 2 >= 0
restricted = []
for s in [F(-1, 2), F(0), F(1, 2)]:
    assert value([s, -s]) == 1 - s * s
    restricted.append({'s': text(s), 'value_on_x_s_y_minus_s': text(value([s, -s]))})

# J(theta)=min(theta,0); its negative stratum has f(theta)=theta, K=0.
# On C=[-h,h], c=0, the correct stratum-restricted linear supremum is0.
# Taking the linear supremum on all C instead yields h, a first-order error.
linear_rows = []
for n in range(1, 7):
    h = F(1, 1 << n)
    correct, enlarged, true_max, K = F(0), h, F(0), F(0)
    assert correct == true_max and enlarged - true_max == h and K == 0
    linear_rows.append({'h': text(h), 'K': '0', 'true_max': '0',
                        'stratum_linear_upper': '0',
                        'full_cell_linear_upper': text(enlarged),
                        'extra_error_over_h': '1'})

result = {
    'status': 'passed', 'scope': 'exact analytic toy identities only; no candidate halo computation',
    'translation': {'periods': ['10', '10'], 'lambda': '1', 'mass': '4',
                    'center': ['5+s', '5'], 'R': '2', 't': '1',
                    'uniform_tail_bound': text(tail),
                    'valid_parameter_range': '|s|<=1/2',
                    'point_area_for_every_s': '4', 'point_ratio_for_every_s': '1',
                    'continuous_family_union_area': '4+4h',
                    'family_union_ratio_minus_true_max': 'h', 'rows': rows},
    'bilinear': {'function': '1+xy', 'H': [['0', '1'], ['1', '0']],
                 'exact_eigenpairs': [{'vector': ['1', '1'], 'eigenvalue': '1'},
                                      {'vector': ['1', '-1'], 'eigenvalue': '-1'}],
                 'rational_unit_positive_direction': list(map(text, positive)),
                 'rational_unit_positive_gain': text(value(positive) - 1),
                 'rational_unit_negative_direction': list(map(text, negative)),
                 'rational_unit_negative_gain': text(value(negative) - 1),
                 'unit_ball_max_gain': '1/2',
                 'max_gain_proof': '2xy<=x^2+y^2<=1; equality at x=y with x^2=y^2=1/2',
                 'restricted_path_identity': 'J(s,-s)=1-s^2',
                 'restricted_path_rows': restricted},
    'stratum_linear_scope': {'function': 'min(theta,0)',
                             'negative_stratum_extension': 'f(theta)=theta',
                             'C': '[-h,h]', 'c': '0', 'K': '0',
                             'rows': linear_rows},
    'evidence_boundary': 'Formula identities establish the toy claims; tabulated samples alone are not used as universal proofs.',
    'python': sys.version, 'platform': platform.platform(),
    'cpu_seconds': time.process_time() - start,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'actual_command': [sys.executable, str(Path(__file__).resolve()),
                       str(Path(sys.argv[1]).resolve())]}
Path(sys.argv[1]).write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items()
                  if k not in ['translation', 'bilinear', 'stratum_linear_scope']}, indent=2))
