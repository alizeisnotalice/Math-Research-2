"""Bracket the area of the R105 AH cap with chords and pairs of tangent lines."""
from fractions import Fraction
import json
from pathlib import Path
import sys

from flint import arb

source = Path(__file__).with_name('certify_ah_eh_combined.py')
prefix, _ = source.read_text().split('assert cut(q(wf))>0 and cut(z)<0')
geometry = {'__file__': str(source)}
argv = sys.argv
try:
    sys.argv = [str(source), '--combined']
    exec(compile(prefix, str(source), 'exec'), geometry)
finally:
    sys.argv = argv

position = geometry['q']
cut = geometry['cut']
scalar = geometry['A']
add, sub, mul, cross = (geometry[k] for k in ('add', 'sub', 'mul', 'cross'))
pivot = geometry['pivot']


def curve_with_tangent(parameter):
    points, tangents = geometry['point_jets'](parameter)
    return geometry['reflection'](points[2]), geometry['reflection'](tangents[2])


lo, hi = Fraction(37, 1000), Fraction(38, 1000)
assert cut(position(lo)) < 0 < cut(position(hi))
for _ in range(80):
    mid = (lo + hi) / 2
    value = cut(position(mid))
    if value < 0:
        lo = mid
    else:
        assert value > 0
        hi = mid
root = geometry['box'](lo, hi)
parts = 512
left, left_tangent = curve_with_tangent(scalar(0))
lower, extra = arb(0), arb(0)
for i in range(1, parts + 1):
    right, right_tangent = curve_with_tangent(root * i / parts)
    chord = sub(right, left)
    denominator = cross(left_tangent, right_tangent)
    assert denominator > 0
    forward = cross(chord, right_tangent) / denominator
    backward = cross(left_tangent, chord) / denominator
    assert forward > 0 and backward > 0
    tip = add(left, mul(forward, left_tangent))
    triangle = cross(sub(tip, left), chord) / 2
    assert triangle > 0
    lower += cross(sub(left, pivot), sub(right, pivot)) / 2
    extra += triangle
    left, left_tangent = right, right_tangent
upper = lower + extra
green_file = Path(__file__).with_name('ah-eh-combined-area-360.json')
green = json.loads(green_file.read_bytes())['AH_cap_area_exact_interval']
a, b = map(lambda t: scalar(Fraction(t)), green)
assert 0 < lower < a and b < upper and extra < scalar(Fraction(1, 10**13))
print(json.dumps(dict(parts=parts, chord_lower=str(lower), tangent_upper=str(upper),
                      bracket_width=str(extra), green_inside=True), indent=2))
