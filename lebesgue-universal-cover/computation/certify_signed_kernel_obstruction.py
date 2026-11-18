"""Bound the full balanced capacity of each explicit common inner polygon."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import sys

from flint import arb, ctx

POSE_SCALE = 2**44
FIT_SCALE = 10**12
TARGET = Fraction('.83743')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def round_up(value):
    return Fraction((value*FIT_SCALE).__ceil__(), FIT_SCALE)


def run(comp, scan_path, core_path, witnesses_path):
    sys.path.insert(0, str(comp.resolve()))
    from audit_mixed_support_blind_box import support_upper, vertices, q, dot
    from certify_five_harmonic import placement_area
    ctx.prec = 192
    scan = load(scan_path)
    core = load(core_path)
    witnesses = load(witnesses_path)
    assert core['source_scan_sha256'] == digest(scan_path)
    assert scan['key'] == core['tree_key'] == [1202, 32259]
    assert scan['classification_axis'] == core['classification_axis'] == 3
    assert scan['target'] == core['target'] == str(TARGET)
    points_path = comp/'finite-point-leaf-input.json'
    assert digest(points_path) == 'b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc'
    data = load(points_path)
    points = [tuple(Fraction(a, data['denominator']) for a in p) for p in data['points']]
    assert set(points) == {(x, -y) for x, y in points}
    assert all((x-u)**2+(y-v)**2 < 1 for x, y in points for u, v in points)
    output = []
    for index, (model, boxes, witness) in enumerate(zip(
        core['common_inner_certificates'], scan['live_boxes'], witnesses, strict=True
    )):
        ref = witness['ref']
        assert all(type(x) is int for x in ref)
        assert any(all(lo <= x <= hi for lo, x, hi in zip(
            box[:11], ref, box[11:], strict=True
        )) for box in boxes)
        assert scan['boundaries'][index] <= ref[3] <= scan['boundaries'][index+1]
        pose = [Fraction(a, POSE_SCALE) for a in ref]
        angle = Fraction(witness['angle'], FIT_SCALE)
        offsets = [Fraction(x, FIT_SCALE) for x in witness['translations']]
        assert len(offsets) == 12
        si, co = q(angle).sin_cos()
        placed_points = [
            (q(offsets[-2])+co*q(x)-si*q(y), q(offsets[-1])+si*q(x)+co*q(y))
            for x, y in points
        ]
        shapes = [
            vertices(n, q(pose[j-1]) if j else arb(0))
            for j, n in zip((0, 3, 6, 9), (3, 5, 7, 9), strict=True)
        ]
        poly = [tuple(map(Fraction, p)) for p in model['polygon']]
        mixed = Fraction(0)
        for u, v in zip(poly, poly[1:]+poly[:1], strict=True):
            n = (v[1]-u[1], u[0]-v[0])
            length = q(dot(n, n)).sqrt()
            normal = tuple(map(q, n))
            shifted_disk = q(dot(offsets[:2], n)) + length/2
            support = max(
                q(dot(u, n)).upper(),
                shifted_disk.upper(),
                max(dot(p, normal).upper() for p in placed_points),
                *((q(dot(offsets[j:j+2], n))+support_upper(shape, normal, length)).upper()
                  for j, shape in zip((2, 4, 6, 8), shapes, strict=True)),
            )
            mixed += Fraction(int((support*10**15).ceil().unique_fmpz()), 2*10**15)
        a5 = Fraction(placement_area(pose, arb.pi(), k=1024)['upper'])
        bound = round_up(mixed)
        assert bound < TARGET
        if index:
            assert a5 < TARGET
        output.append(dict(
            class_number=index,
            five_area_upper=str(round_up(a5)),
            genuine_five_sublevel=a5 < TARGET,
            mixed_budget_max=str(bound),
            strict_gap=str(TARGET-bound),
        ))
    return dict(
        classes=output,
        target=str(TARGET),
        source_sha256=digest(Path(__file__)),
        witnesses_sha256=digest(witnesses_path),
        p30_sha256=digest(points_path),
        signed_scan_sha256=digest(scan_path),
        shared_polygons_sha256=digest(core_path),
    )


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions required.')
    parser = ArgumentParser(description=__doc__)
    for name in ('comp', 'scan', 'cores', 'witnesses'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.comp, args.scan, args.cores, args.witnesses),
                     separators=(',', ':')))
