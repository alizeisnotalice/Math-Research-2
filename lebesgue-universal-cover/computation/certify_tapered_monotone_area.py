"""Area forced by a thirteen-arc input for supersets of the R203 host, with adjustable strip slack."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from flint import arb, ctx
from audit_eh_continuous_placement import project
from certify_ah_tapered_minkowski_closure import geometry


def run(comp):
    env, geometry_sha = geometry(comp)
    ctx.prec = 160
    g = env['g']
    n, m, c3, b3, *_ = g.setup()
    e3 = g.line_pair(n[4], m[3])
    f3 = g.line_pair(n[3], m[2])
    # Endpoints are the geometric intersections/projections defining the two Gibbs
    # tilted faces and the transverse endpoints of the new deep cut (R203).
    segments = {
        'I-E3': (project(b3, m[3]), e3),
        'K-C3': (project(f3, m[5]), c3),
        'deep': (env['W'], env['crosspoint']),
    }
    wall = {name:(normal,height) for name,normal,height in env['active']}
    lengths = {}
    margin = arb(1)
    for name, (a,b) in segments.items():
        u,h = wall[name]
        assert all((h-g.dot(u,p)).contains(0) for p in (a,b))
        direction = g.minus(b,a)
        cut = arb(1)/10**6
        for t in (cut, 1-cut):
            pt = g.plus(a,g.scale(t,direction))
            for other,v,offset in env['active']:
                if name != other:
                    gap = offset-g.dot(v,pt)
                    assert gap > 0, (name,other)
                    margin = min(margin,gap.lower())
        lengths[name] = ((1-2*cut)*g.norm(direction)).lower()
    mirror_file = comp/'ah-tapered-mirror-counterexample-audit.json'
    profile = json.loads(mirror_file.read_bytes())
    arc = profile['binary_curvature_Reuleaux_polygon']
    assert arc['same_compatibility_roots'] and arc['full_period_unit_arc_count'] == 13
    assert [x['rotation_degrees'] for x in arc['six_rotations']] == list(range(0,360,60))
    cells = profile['compatibility_closed_cells_half_chamber']
    assert cells[0]['domain_degrees'] == ['0','15/2'] and cells[-1]['domain_degrees'][1]=='30'
    speed = -arb(cells[0]['upper'])
    drop = min(-arb(c['upper']) for c in cells[1:])-3*arb(arc['Hausdorff_error_upper'])
    assert speed > arb('.11248') and drop > arb('.015')
    area = env['area']
    threshold = arb('0.844082046521218581802948589770767704354')
    r = 1/arb(3).sqrt()
    rows = []
    for slack in ('0','.0000002'):
        eta = arb(slack)
        assert 3*eta < drop
        angle = 3*eta/speed
        motion = r*angle+2*r*(eta+r*angle)
        costs = []
        for pose in arc['six_rotations']:
            row = pose['true_face']
            clearance = arb(pose['violation_lower']) - motion
            assert clearance>0
            costs.append((lengths[row]*clearance/2).lower())
        bound = (area+min(costs)).lower()
        assert bound > threshold
        rows.append(dict(per_side_widening=slack,angle_limit=str(angle.upper()),
                         motion_upper=str(motion.upper()),triangle_areas_lower=list(map(str,costs)),
                         total_area_lower=str(bound),excess_over_B=str((bound-threshold).lower())))
    return dict(geometry_sha256=geometry_sha,mirror_audit_sha256=sha256(mirror_file.read_bytes()).hexdigest(),
                segment_lengths_lower={label:str(x) for label,x in lengths.items()},
                other_wall_minimum=str(margin),cases=rows)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions must be enabled.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('computation',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    result=run(args.computation.resolve())
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for row in result['cases']:
        print(row['per_side_widening'],row['total_area_lower'],row['excess_over_B'])
