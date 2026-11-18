#!/usr/bin/env python3
"""Exclude the widened-strip host using all completions of one 18-point witness."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_general_relocation import exact_fraction, rotate_angle
from certify_eh_no_interior_roots import verify
from certify_eh_polyline_obstruction import witness_points


STRIP_LABELS = ('C3-B1', 'B1-Q', 'P-F1', 'F1-H', 'E3-D1', 'D1-G')
DILATION = Fraction(1, 2000000)


def verify_host(host, case, dilation=DILATION):
    result = verify(case)
    points, _ = witness_points(case)
    normals, *_ = g.setup()
    den = host['planes_denominator']
    by_name = {name:(a,b,c) for name,a,b,c in host['rational_planes']}
    depths = list(map(Fraction,host['outward_strip_allowances_exact']))
    assert len(depths) == 6 and min(depths) >= 0
    for i,name in enumerate(STRIP_LABELS):
        a,b,c = by_name[name]
        rounded = [exact_fraction(Fraction(x,den)) for x in (a,b)]
        error = sum(abs(x-y) for x,y in zip(rounded,normals[i]))
        assert exact_fraction(Fraction(c,den))+error < g.HALF+exact_fraction(depths[i])

    # For C+B(dilation) each planar support is allowed dilation more.
    cap = max(sum(depths[i] for i in (0,2,4)),sum(depths[i] for i in (1,3,5)))+3*dilation
    interior = g.arb(result['uniform_F_lower_interior'])
    slope = g.arb(result['uniform_F_prime_lower_first_quarter_degree'])
    assert interior > exact_fraction(cap) and slope > 0
    angle = exact_fraction(cap)/slope
    assert angle < g.PI/720

    # The old witness fixes all original hexagon supports at 1/2. Estimate
    # motion from the nearest 60-degree position using the n1,n2 strips.
    radius = 1/g.arb(3).sqrt()
    pair_strip = exact_fraction(max(depths[i] for i in (1,2,4,5))+dilation)
    displacement = 3*radius*angle+2*pair_strip

    edges = []
    for name,a,b,c in host['rational_planes']:
        normal = g.arb(a),g.arb(b)
        length = g.norm(normal)
        edges.append((name,g.scale(1/length,normal),g.arb(c)/length))
    gaps = []
    for j in range(6):
        for reflect in (False,True):
            corner_angle = j*g.PI/3
            image = {name:rotate_angle((x,-y if reflect else y),corner_angle) for name,(x,y) in points.items()}
            candidates = ((g.dot(u,p)-b,edge,name) for edge,u,b in edges for name,p in image.items())
            violation, edge, name = max(candidates,key=lambda row:row[0].lower())
            gap = violation-displacement-exact_fraction(dilation)
            assert gap > 0, (j,reflect,edge,name,str(gap))
            gaps.append((j,reflect,edge,name,violation.lower(),gap.lower()))

    # For every completion choose support points at 5000 equal angles.
    # The support-point Lipschitz estimate makes their hull an inner net.
    n = 5000
    assert (g.PI/n)**2 + g.arb('0.00000000023') < exact_fraction(dilation)
    return {'degree_free_orientation_argument':'even allowing the stated outward Euclidean neighborhood',
            'dilation_exact':str(dilation),'original_witness':result,
            'required_sector_compatibility_abs_upper_exact':str(cap),
            'nearest_hex_symmetry_angle_radians_upper':str(angle.upper()),
            'point_displacement_upper':str(displacement.upper()),
            'least_overflow_after_displacement_and_dilation_lower':str(min(row[-1] for row in gaps)),
            'hex_symmetry_violations':[list(row[:4])+[str(row[4]),str(row[5])] for row in gaps],
            'existence_of_rational_finite_obstruction_at_most_n_points':n}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('host',type=Path)
    parser.add_argument('witness',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    host, case = (json.loads(path.read_bytes()) for path in (args.host,args.witness))
    result = verify_host(host,case)
    result['input_sha256'] = {path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in (args.host,args.witness)}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
