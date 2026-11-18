#!/usr/bin/env python3
"""Rigorous width-one witness tests for two-segment EH boundaries."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_all_completions import Lens
from certify_eh_general_relocation import exact_fraction
from certify_eh_no_interior_roots import verify


def clip_geometry(params):
    first, second, height = params
    assert Fraction(first) > Fraction(second) > 0 and 0 <= Fraction(height) < 1
    enclosing = g.trial(first, tight=True)
    g.trial(second, tight=True)
    _, _, center, lower_center, _, upper, _, lower, _ = g.setup()

    def move(depth):
        cosine, sine = g.unit(exact_fraction(Fraction(depth)))
        radial = g.minus(upper, center)
        return g.plus(center, g.minus(g.scale(cosine, radial), g.scale(sine, g.rot(radial))))

    first_upper, second_upper = move(first), move(second)
    hinge = g.plus(second_upper, g.scale(exact_fraction(Fraction(height)), g.minus(lower, second_upper)))
    edges = g.minus(hinge, first_upper), g.minus(lower, hinge)
    turn = g.cross(*edges)
    assert all(v[1] < 0 for v in edges) and turn > 0
    assert g.cross(g.rot(g.minus(first_upper, center)), edges[0]) > 0
    assert g.cross(edges[1], g.rot(g.minus(lower, lower_center))) > 0
    normals = [g.scale(1 / g.norm(v), g.neg(g.rot(v))) for v in edges]
    offsets = [g.dot(normal, point) for normal, point in zip(normals, (first_upper, hinge))]
    assert g.dot(normals[0], lower) < offsets[0] and g.dot(normals[1], first_upper) < offsets[1]
    gain = g.arb(enclosing['gain_area']) - turn / 2
    return normals, offsets, gain, first_upper, hinge


def witness_points(case):
    n, m, c3, *_ = g.setup()
    phi = exact_fraction(Fraction(case['support_angle_degrees'])) * g.PI / 180
    axis = n[5]

    def reflect(point):
        return g.minus(g.scale(2 * g.dot(axis, point), axis), point)

    raw = list(map(g.arb, case['coordinates']))
    z, zq = tuple(raw[:2]), tuple(raw[2:4])
    b, bq = (g.plus(c3, (g.arb(0), height)) for height in raw[4:6])
    w, wq = tuple(raw[6:8]), tuple(raw[8:10])
    q = g.plus(g.scale(g.HALF, axis), g.scale(raw[10], g.rot(axis)))
    points = {'z': z, 'b': b, 'b_op': g.minus(b, n[0]), 'w': w, 'w_op': g.plus(w, m[5]),
              'q': q, 'q_op': g.minus(q, axis), 'z_Q': reflect(zq), 'b_Q': reflect(bq),
              'b_Q_op': reflect(g.minus(bq, n[0])), 'w_Q': reflect(wq), 'w_Q_op': reflect(g.plus(wq, m[5]))}
    for j, i in enumerate((11, 13, 15)):
        point = tuple(raw[i:i+2])
        direction = g.scale((-1)**j, g.unit(phi + j*g.PI/3))
        points[f'contact_{j}'] = point
        points[f'contact_{j}_op'] = g.minus(point, direction)
    return points, reflect


def test_shape(witness, params):
    normals, offsets, gain, _, _ = clip_geometry(params)
    points, reflect = witness_points(witness)
    lens = Lens(points)
    margins = []
    for enabled in (False, True):
        directions = list(map(reflect, normals)) if enabled else normals
        violations = [[(g.dot(normal, point)-offset, label, i) for label, point in points.items()]
                      for i, (normal, offset) in enumerate(zip(directions, offsets))]
        v, label, edge = max((item for row in violations for item in row), key=lambda r: r[0].lower())
        margins.append({'pose': 'Q2' if enabled else 'identity',
                        'largest_violation_lower': str(v.lower()), 'witness_point': label, 'hinge_edge': edge,
                        'outer_lens_slack_lower': [str((b-lens.upper_support(u)).lower())
                                                  for u, b in zip(directions, offsets)],
                        'finite_set_excluded': bool(v > 0)})
    cprime = g.unit(11*g.PI/6-g.SIGMA)
    old_ah_clearance = g.HALF-lens.upper_support(reflect(cprime))
    return {'parameters': params, 'EH_gain_lower': str(gain.lower()), 'EH_gain_upper': str(gain.upper()),
            'poses': margins, 'reflected_Cprime_slack_lower': str(old_ah_clearance.lower())}


def main(args):
    data = json.loads(args.input.read_bytes())
    fixed = data['reference_witness']
    baseline = g.arb(data['baseline_gain_upper'])
    fixed_certificate = verify(fixed)
    reference_shape = test_shape(fixed, data['restored_reference_shape'])
    assert g.arb(reference_shape['EH_gain_lower']) > baseline
    assert all(g.arb(s) > 0 for s in reference_shape['poses'][1]['outer_lens_slack_lower'])
    assert g.arb(reference_shape['reflected_Cprime_slack_lower']) > 0

    ruled_out = []
    for params in data['screened_better_shapes']:
        row = test_shape(fixed, params)
        assert g.arb(row['EH_gain_lower']) > baseline
        if all(pose['finite_set_excluded'] for pose in row['poses']):
            ruled_out.append({'obstruction': 'reference', **row})

    for obstruction in data['new_witnesses']:
        params = obstruction['polyline']
        assert Fraction(params[0]) == Fraction(obstruction['delta'])
        proof = verify(obstruction)
        result = test_shape(obstruction, params)
        assert all(pose['finite_set_excluded'] for pose in result['poses'])
        ruled_out.append({'obstruction': 'additional', **result,
                          'F_inner_sector_lower': proof['uniform_F_lower_interior'],
                          'F_inner_cells': proof['interior_angle_cell_count'],
                          'squared_diameter_slack': proof['squared_slack_other_pairs']})

    checked = {tuple(map(Fraction, row['parameters'])) for row in ruled_out}
    assert checked >= {tuple(map(Fraction, params)) for params in data['screened_better_shapes']}

    baseline_lower = g.arb(data['baseline_gain_lower'])
    angular_walls = []
    previous = Fraction(0)
    for case in data['hinge_interval_witnesses']:
        params, interval = case['polyline'], case['wall_hinge_range']
        first, second, hinge = params
        assert Fraction(second) == Fraction('0.00055') and Fraction(first) == Fraction(case['delta'])
        assert Fraction(hinge) == Fraction(interval[0]) == previous < Fraction(interval[1]) < 1
        proof = verify(case)
        normals, offsets, _, upper, kink = clip_geometry(params)
        values = list(map(g.arb, case['coordinates']))
        witnesses = tuple(values[:2]), tuple(values[2:4])
        distances = []
        for point, part in zip(witnesses, case['branch']):
            assert part in (0, 1)
            if part == 0:
                assert kink[1] < point[1] < upper[1]
            distance = g.dot(normals[part], point)-offsets[part]
            assert distance > 0
            distances.append(str(distance.lower()))
        _, _, area, _, _ = clip_geometry([first, second, interval[1]])
        assert area < baseline_lower
        angular_walls.append({'hinge_fraction': interval, 'critical_upper_depth': first,
                              'area_gain_upper_if_below_critical_depth': str(area.upper()),
                              'F_inner_sector_lower': proof['uniform_F_lower_interior'],
                              'both_ear_violations_lower': distances})
        previous = Fraction(interval[1])
    assert previous == Fraction('0.99')
    assert Fraction(data['hinge_interval_witnesses'][-1]['delta']) == Fraction('0.00055222')
    _, _, center, _, _, old_upper, *_ = g.setup()
    old_radial = g.minus(old_upper, center)
    cap = exact_fraction(Fraction('0.006'))
    assert old_radial[0] < 0 and -g.arb('0.008') < old_radial[1] < -g.arb('0.007')
    assert cap.cos() > g.arb('0.9999') and cap.sin() < g.arb('0.0061')
    last_single_gain = g.arb(g.trial('0.00055222', tight=True)['gain_area'])
    remaining_gain = last_single_gain-baseline_lower
    assert 0 < remaining_gain < g.arb('4.52e-12')

    output = {'reference_F_inner_sector_lower': fixed_certificate['uniform_F_lower_interior'],
              'reference_witness_class_fits': reference_shape,
              'strictly_improved_shapes_ruled_out': len(checked), 'ruled_out': ruled_out,
              'hinge_intervals_at_fixed_second_depth': angular_walls,
              'remaining_possible_area_gain_over_single_chord_upper': str(remaining_gain.upper())}
    args.output.write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    main(parser.parse_args())
