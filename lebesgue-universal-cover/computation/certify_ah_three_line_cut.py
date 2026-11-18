#!/usr/bin/env python3
"""Reuse the two-line AH proof for an additional cut at its interior kink."""
import argparse
import contextlib
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

from flint import arb, ctx


SOURCE_SHA256 = '3d3781ceccaf14b9db4c0bb096e6d24a98c825c0a5f1eadbe6df639a57c86e9b'
OLD_CUTS = 'cuts=[(F("0.355"), F("0.00060")), (F("0.410"), F("0.00030"))]'


def specialized_source(path, pivot, angle):
    raw = path.read_bytes()
    assert sha256(raw).hexdigest() == SOURCE_SHA256
    text = raw.decode()
    extra = f', (F("{pivot}"), F("{angle}"))]'
    replacements = [(OLD_CUTS,OLD_CUTS[:-1]+extra),
                    ('eh_depth = "0.00055"','eh_depth = "0.000555"'),
                    ('"eh-cap-switch-joint-550.json"','"eh-two-line-cover-certificate.json"'),
                    ('baseline["validated_duals"] == 60',
                     'baseline["validated_duals"] == 63 and baseline["delta_is_enclosing_placement_depth"]')]
    for old,new in replacements:
        assert text.count(old)==1
        text = text.replace(old,new)
    return text


def dyadic(ball):
    mantissa,exponent = map(int,ball.man_exp())
    return Fraction(mantissa)*Fraction(2)**exponent


def format_area(lower, upper, digits=39):
    scale = 10**digits
    a,b = lower*scale//1, -(-upper*scale//1)
    return [f'{i//scale}.{i%scale:0{digits}d}' for i in (a,b)]


def triangle(ns):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_polyline_obstruction import clip_geometry

    z,D,n,length = (ns[k] for k in ('z','D','normal','length'))
    (p1,_),(p2,_),(p3,_) = ns['cuts']
    b1,b2,b3 = ns['tangents']
    ar = ns['A']
    a1,a2,a3 = (ar(p)*b for p,b in zip((p1,p2,p3),(b1,b2,b3)))
    assert b1>b3>b2>0
    w0 = (a1-a2)/(b1-b2)
    w1 = (a1-a3)/(b1-b3)
    w2 = (a3-a2)/(b3-b2)
    assert 0 < w1 < w0 < w2 < ar(p2)
    excess = a3-b3*w0-a1+b1*w0
    assert excess > 0

    def point(w, height):
        return g.plus(g.plus(z,g.scale(w,D)),g.scale(length*height,n))

    corners = [point(w,a-b*w) for w,a,b in ((w1,a1,b1),(w0,a1,b1),(w2,a2,b2))]
    g.SIGMA = ar(Fraction('1.5497'))*g.PI/180
    inner,_ = polygon('0.0005522',32)
    directions, offsets, *_ = clip_geometry(['0.000555','0.0005522','0.0075'])
    inner += [(f'R126_EH_{i}',u,b) for i,(u,b) in enumerate(zip(directions,offsets))]
    clear = min((b-g.dot(u,x)).lower() for _,u,b in inner for x in corners)
    assert clear > 0
    area = length**2*(w2-w1)*excess/2
    return {'left_intersection_w':str(w1),'old_kink_w':str(w0),'right_intersection_w':str(w2),
            'triangle_strict_old_curve_clearance_lower':str(clear),
            'area_gain_exact':[str(dyadic(v)) for v in (area.lower(),area.upper())]},area,inner


def fit_archived_witnesses(archive, ns, inner):
    import audit_optimized_eh_chord as g

    rows = list(inner)
    for i,(direction,pivot) in enumerate(zip(ns['cut_dirs'],ns['pivot_set'])):
        normal = g.neg(g.rot(direction))
        rows.append((f'AH_{i}',normal,g.dot(normal,pivot)))
    verified = {}
    for name,count in (('widened-projected-obstruction-points.json',60),
                       ('widened-projected-obstruction-18-points.json',18)):
        raw = (archive/name).read_bytes()
        data = json.loads(raw)
        den = data['denominator']
        points = [tuple(ns['A'](Fraction(x,den)) for x in p) for p in data['points']]
        assert len(points)==count
        values = [((b-g.dot(u,p)).lower(),label,j) for label,u,b in rows for j,p in enumerate(points)]
        clearance,edge,index = min(values,key=lambda row:row[0])
        assert clearance > g.arb('2.15e-14')
        verified[name] = {'all_points_fit_in_identity_pose':True,'active_edge':edge,'active_point_index':index,
                          'clearance_lower':str(clearance),'input_sha256':sha256(raw).hexdigest()}
    return verified


def separate_widened_host(archive, ns):
    import audit_optimized_eh_chord as g
    from certify_eh_polyline_obstruction import clip_geometry

    _,_,_,upper,_ = clip_geometry(['0.000555','0.0005522','0.0075'])
    for along,pivot in zip(ns['cut_dirs'],ns['pivot_set']):
        assert g.cross(along,g.minus(upper,pivot)) > 0
    path = archive/'r112-widened-host-certificate.json'
    failed_host = json.loads(path.read_bytes())
    a,b,c = next((a,b,c) for name,a,b,c in failed_host['rational_planes'] if name=='U1-Te')
    direction = g.arb(a),g.arb(b)
    gap = (g.dot(direction,upper)-g.arb(c))/g.norm(direction)
    assert gap > g.arb('0.00021')
    return {'retained_R126_EH_upper_end':list(map(str,upper)),
            'gap_from_depth_0p0072_chord_lower':str(gap.lower()),
            'excluded_host_sha256':sha256(path.read_bytes()).hexdigest()}


def run(archive, pivot, angle):
    sys.path.insert(0,str(archive))
    old_script = archive/'certify_ah_two_line_cap.py'
    text = specialized_source(old_script,pivot,angle)
    ns={'__file__':str(old_script),'__name__':'__main__'}
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        exec(compile(text,str(old_script),'exec'),ns)
    placement = json.loads(captured.getvalue())
    assert placement['passed'] and placement['EH_depth_radians_exact']=='0.000555'
    assert all(arb(placement[k]) > 0 for k in ('uniform_perturbed_compatibility_slack_lower',
                                               'uniform_actual_placement_slack_lower'))

    ctx.prec = 160
    shape, area, inner = triangle(ns)
    prior_path = archive/'ah-eh-two-line-combined.json'
    eh_path = archive/'eh-two-line-cover-certificate.json'
    prior = json.loads(prior_path.read_bytes())
    assert prior['passed'] and prior['new_EH_depth_exact']=='0.000555'
    assert prior['input_sha256'][eh_path.name]==sha256(eh_path.read_bytes()).hexdigest()
    a,b = map(Fraction,prior['combined_candidate_area_exact_interval'])
    lo,hi = a-dyadic(area.upper()), b-dyadic(area.lower())
    assert lo < hi < a

    return {'new_AH_third_pivot_exact':str(pivot),'new_AH_third_angle_exact':str(angle),
            'three_line_AH_placement':placement,'triangle':shape,
            'separation_from_R145_host':separate_widened_host(archive,ns),
            'witnesses_at_identity':fit_archived_witnesses(archive,ns,inner),
            'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'derived_AH_proof_sha256':sha256(text.encode()).hexdigest(),
            'baseline_area_sha256':sha256(prior_path.read_bytes()).hexdigest()}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('pivot',type=Fraction)
    parser.add_argument('angle',type=Fraction)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    assert 0 < args.pivot < Fraction('.410') and 0 < args.angle < Fraction('.00060')
    certificate = run(args.archive.resolve(),args.pivot,args.angle)
    args.output.write_text(json.dumps(certificate,indent=2)+'\n')
