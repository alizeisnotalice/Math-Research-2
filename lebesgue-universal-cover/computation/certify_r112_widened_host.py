#!/usr/bin/env python3
"""Fit all completions of the R112 failure witness into a widened-strip host."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from audit_eh_continuous_placement import polygon
from certify_eh_all_completions import Lens
from certify_eh_general_relocation import exact_fraction, rotate_angle
from certify_eh_no_interior_roots import verify
from certify_eh_polyline_obstruction import witness_points


WIDTH_LABELS = ('C3-B1', 'B1-Q', 'P-F1', 'F1-H', 'E3-D1', 'D1-G')
HAT_LABELS = ('I-E3', 'K-C3')
DENOMINATOR = 10**20
SUBDIVISIONS = 40
DEPTH = '0.0072'
PAD = Fraction('0.000000001')
EAST_EXTENSION = Fraction('0.00000130')
HAT_ROUNDING_PAD = Fraction('0.000000000000000001')
TRANSLATION = tuple(map(Fraction, ('0.003358657267784', '0.011268447942153')))


def dyadic(ball):
    mantissa, exponent = map(int, ball.man_exp())
    return Fraction(mantissa) * Fraction(2)**exponent


def ah_lines():
    n, m, *_ = g.setup()
    e3 = g.line_pair(n[4], m[3])
    c = -e3[0], e3[1]
    root3 = g.arb(3).sqrt()
    mid = -g.HALF/2, -root3/4
    l0x = (root3-1)/2
    z = g.circle_pair((l0x, (l0x-1)/root3), mid, -1)
    v = g.minus(z, c)
    squared = g.dot(v, v)
    ta = g.plus(c, g.scale(1/squared, g.minus(v, g.scale((squared-1).sqrt(), g.rot(v)))))
    direction = g.minus(ta, z)
    along = g.scale(-1/g.norm(direction), direction)
    lines = []
    for i, (p, a) in enumerate((('0.355','0.00060'), ('0.410','0.00030'))):
        pivot = g.plus(z,g.scale(g.arb(p),direction))
        normal = g.neg(g.rot(rotate_angle(along,g.arb(a))))
        lines.append((f'AH_{i}', normal, g.dot(normal,pivot)))
    return lines


def rational_planes():
    unit_width, _ = polygon(DEPTH, SUBDIVISIONS)
    planes = []
    for label, normal, bound in unit_width + ah_lines():
        if label in WIDTH_LABELS:
            bound += exact_fraction(EAST_EXTENSION if label == WIDTH_LABELS[0] else PAD)
        elif label in HAT_LABELS:
            # The two circle-tangent inequalities need a rounding allowance.
            bound += exact_fraction(HAT_ROUNDING_PAD)
        a = tuple((dyadic(v.mid())*DENOMINATOR)//1 for v in normal)
        b = (dyadic(bound.lower())*DENOMINATOR)//1-3
        # On the starting square, each rational plane lies inside the exact
        # geometric plane even after accounting for direction rounding.
        error = sum(abs(v-exact_fraction(Fraction(w,DENOMINATOR))) for v,w in zip(normal,a))
        assert exact_fraction(Fraction(b,DENOMINATOR))+error < bound
        planes.append((label, *a, b))
    return planes


def polytope_area(planes):
    vertices = [(Fraction(x),Fraction(y)) for x,y in ((-1,-1),(1,-1),(1,1),(-1,1))]
    for _, a, b, limit in planes:
        new = []
        for x, y in zip(vertices, vertices[1:]+vertices[:1]):
            delta = a*x[0]+b*x[1]-limit
            other = a*y[0]+b*y[1]-limit
            if delta <= 0:
                new.append(x)
            if delta*other < 0:
                ratio = delta/(delta-other)
                new.append(tuple(u+(v-u)*ratio for u,v in zip(x,y)))
        vertices = new
        assert len(vertices) >= 3
    return sum(x[0]*y[1]-x[1]*y[0] for x,y in zip(vertices, vertices[1:]+vertices[:1]))/2, len(vertices)


def certify_placement(planes, lens, angle, translation):
    minimum = {}
    square = [(f'bounding_square_{i}',a*DENOMINATOR,b*DENOMINATOR,DENOMINATOR)
              for i,(a,b) in enumerate(((1,0),(-1,0),(0,1),(0,-1)))]
    for label, a, b, c in planes + square:
        normal = tuple(exact_fraction(Fraction(x,DENOMINATOR)) for x in (a,b))
        length = g.norm(normal)
        u = g.scale(1/length,normal)
        clearance = exact_fraction(Fraction(c,DENOMINATOR))/length - g.dot(u,translation)
        clearance -= lens.upper_support(rotate_angle(u,-angle)) if lens else g.HALF
        assert clearance > 0, (label, str(clearance))
        minimum[label] = clearance.lower()
    return minimum


def decimal_bounds(value, digits):
    scale = 10**digits
    lower, upper = value*scale//1, -(-value*scale//1)
    return [f'{i//scale}.{i%scale:0{digits}d}' for i in (lower,upper)]


def run(input_path, baseline_path, ah_area_path, comparison_path):
    cases = json.loads(input_path.read_bytes())
    case = next(c for c in cases if Fraction(c['delta']) == Fraction('0.0005525'))
    assert Fraction(case['sigma']) == Fraction('1.5497')
    assert Fraction(case['support_angle_degrees']) == 36
    original = verify(case)
    points, _ = witness_points(case)
    reflected_lens = Lens({name:(x,-y) for name,(x,y) in points.items()})
    angle = -84*g.PI/180
    translation = tuple(map(exact_fraction,TRANSLATION))
    planes = rational_planes()
    margins = certify_placement(planes,reflected_lens,angle,translation)
    shape_clearance, worst_shape_edge = min((v,k) for k,v in margins.items() if k not in WIDTH_LABELS)
    width_clearance, worst_width_edge = min((margins[name],name) for name in WIDTH_LABELS)
    assert shape_clearance > g.arb('0.00000140')
    assert width_clearance > g.arb('0.00000000099999')

    center = (g.arb(0),g.arb(0))
    disk = certify_placement(planes,None,g.arb(0),center)
    n, *_ = g.setup()
    triangle = Lens({str(i): g.line_pair(n[i],n[(i+1)%6]) for i in (0,2,4)})
    reuleaux = certify_placement(planes,triangle,g.arb(0),center)

    area, num_vertices = polytope_area(planes)
    enclosure = decimal_bounds(area,36)
    assert Fraction(enclosure[0]) < area < Fraction(enclosure[1])
    source = json.loads(baseline_path.read_bytes())
    ah_area = json.loads(ah_area_path.read_bytes())
    comparison = json.loads(comparison_path.read_bytes())
    assert Fraction(source['sigma_degrees']) == Fraction(ah_area['sigma_degrees_exact']) == Fraction('1.5497')
    assert ah_area['AH_cut_parameters_exact'] == [['71/200','3/5000'], ['41/100','3/10000']]
    area_upper = exact_fraction(Fraction(source['area_interval'][1]))
    area_upper -= exact_fraction(Fraction(ah_area['AH_cap_area_exact_interval'][0]))
    gain_eh = g.arb(g.trial(DEPTH,tight=True)['gain_area'])
    area_upper -= gain_eh

    max_extension = exact_fraction(EAST_EXTENSION)
    strip_allowances = [exact_fraction(EAST_EXTENSION)] + [exact_fraction(PAD)]*5
    added_area = sum((1+4*max_extension-e)*e for e in strip_allowances)/g.arb(3).sqrt()
    added_area += 4*g.arb(2).sqrt()*exact_fraction(HAT_ROUNDING_PAD)
    area_upper += added_area
    assert exact_fraction(area) < area_upper
    baseline_lower = exact_fraction(Fraction(comparison['combined_candidate_area_decimal_interval'][0]))
    assert baseline_lower-area_upper > g.arb('0.00000086')
    return {'old_obstruction': {'depth':case['delta'], 'F_inner_sector_lower':original['uniform_F_lower_interior']},
            'scope':'all constant-width-one completions of the specified R112 point set; also the unit-diameter disk and regular Reuleaux triangle',
            'EH_depth':DEPTH, 'arc_subdivisions':SUBDIVISIONS,
            'outward_strip_allowances_exact':[str(EAST_EXTENSION)]+[str(PAD)]*5,
            'hat_rounding_allowance_exact':str(HAT_ROUNDING_PAD),
            'fixed_witness_placement':{'reflect_across_x_axis':True, 'angle_degrees':'-84',
                                       'translate_exact':list(map(str,TRANSLATION))},
            'min_other_face_clearance_lower':str(shape_clearance), 'critical_other_face':worst_shape_edge,
            'min_hex_clearance_lower':str(width_clearance), 'critical_hex_face':worst_width_edge,
            'disk_clearance_lower':str(min(disk.values())), 'regular_Reuleaux_clearance_lower':str(min(reuleaux.values())),
            'area_decimal_interval':enclosure, 'polygon_vertices':num_vertices,
            'EH_cut_gain_interval':str(gain_eh), 'strip_cost_upper':str(added_area.upper()),
            'area_upper_including_arc_area_budget':str(area_upper.upper()),
            'gain_over_R126_using_area_budget_lower':str((baseline_lower-area_upper).lower()),
            'input_sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in (input_path,baseline_path,ah_area_path,comparison_path)},
            'planes_denominator':DENOMINATOR, 'rational_planes':planes}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('witness',type=Path)
    parser.add_argument('gibbs_area',type=Path)
    parser.add_argument('AH_area',type=Path)
    parser.add_argument('reference_R126_area',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    result = run(args.witness,args.gibbs_area,args.AH_area,args.reference_R126_area)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
