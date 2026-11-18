#!/usr/bin/env python3
"""Transfer one R167 placement-margin sweep to a jointly selected AH envelope."""
import argparse
import contextlib
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

from flint import arb, ctx

from certify_ah_three_line_cut import dyadic, fit_archived_witnesses, format_area, specialized_source


R167_PIVOT = Fraction('0.404')
R167_ANGLE = Fraction('0.00046')
ANCHORED_CUTS = 'cuts=[(F("0.355"), F("0.00060")), (F("0.410"), F("0.00030")), (F("101/250"), F("23/50000"))]'


def replace_line(text, old, new):
    assert text.count(old)==1
    return text.replace(old,new)


def sweep_r167(archive):
    source = archive/'certify_ah_two_line_cap.py'
    code = specialized_source(source,R167_PIVOT,R167_ANGLE)
    tracked = replace_line(code,'counts=0','counts=0\nroot_by_bin={}\nplacement_by_bin={}')
    tracked = replace_line(tracked,'        best_root=min(best_root,*(f.lower() for f in fs))',
                           '        best_root=min(best_root,*(f.lower() for f in fs))\n'
                           '        root_by_bin[j]=min(root_by_bin.get(j,one),*(f.lower() for f in fs))')
    tracked = replace_line(tracked,'            best_actual_edge=min(best_actual_edge,(ms[worst]-gap).lower())',
                           '            best_actual_edge=min(best_actual_edge,(ms[worst]-gap).lower())\n'
                           '            placement_by_bin[j]=min(placement_by_bin.get(j,one),(ms[worst]-gap).lower())')
    ns = {'__file__':str(source),'__name__':'__main__'}
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        exec(compile(tracked,str(source),'exec'),ns)
    assert json.loads(captured.getvalue())['passed']
    return ns, code, tracked


def new_geometry(code, old, extra):
    text = ANCHORED_CUTS[:-1]+''.join(f', (F("{p}"), F("{angle}"))' for p,angle in extra)+']'
    full = replace_line(code,ANCHORED_CUTS,text)
    marker = '\nroot_seeds=['
    assert full.count(marker)==1
    prefix = full.split(marker)[0]
    ns = {'__file__':old['__file__'],'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix,old['__file__'],'exec'),ns)
    return ns, full


def coefficients(ns):
    return [(b,ns['A'](p)*b,i) for i,((p,_),b) in enumerate(zip(ns['cuts'],ns['tangents']))]


def transfer(old, new):
    # The source sweep subtracts 3*1.07 times the chart deficit in each bin.
    ar = new['A']
    L = new['length']
    b_old = sorted(coefficients(old),key=lambda row:row[0].mid(),reverse=True)
    assert all(x[0]>y[0] for x,y in zip(b_old,b_old[1:]))
    excess = arb(0)
    maximum_norm = arb(1)
    for b,a,i in coefficients(new)[len(b_old):]:
        upper,lower = next((p,q) for p,q in zip(b_old,b_old[1:]) if p[0]>b>q[0])
        lam = (b-lower[0])/(upper[0]-lower[0])
        assert 0 < lam < 1
        # Outward unnormalized normal is affine in the slope b.
        correction = L*(a-lam*upper[1]-(1-lam)*lower[1])
        excess = max(excess,correction.upper())
        maximum_norm = max(maximum_norm,(1+b*b).sqrt().upper())

    rows = []
    for i,(was,now) in enumerate(zip(old['deficit_by_bin'],new['deficit_by_bin'],strict=True)):
        if now.is_zero():
            continue
        assert i in old['root_by_bin'] and i in old['placement_by_bin']
        old_root = old['root_by_bin'][i]
        old_pose = old['placement_by_bin'][i]
        error = max(arb(0),(ar('3.21')*(now-was)).upper())
        root = old_root-error
        core_pose = old_pose-error
        added_pose = (core_pose-excess)/maximum_norm
        assert root > 0 and core_pose > 0 and added_pose > 0, (i,str(root),str(core_pose),str(added_pose))
        rows.append({'bin':i,'R167_deficit':str(was),'new_deficit':str(now),
                     'R167_root_margin_lower':str(old_root),'R167_placement_margin_lower':str(old_pose),
                     'extra_deficit_cost_upper':str(error), 'root_margin_lower':str(root.lower()),
                     'placement_margin_lower':str(min(core_pose.lower(),added_pose.lower()))})
    return {'old_AH_parameter_cells':old['counts'],'angle_cells_each':old['angular_cells'],
            'new_line_affine_offset_cost_upper':str(excess),
            'joint_root_slack_lower':str(min(arb(row['root_margin_lower']).lower() for row in rows)),
            'joint_placement_slack_lower':str(min(arb(row['placement_margin_lower']).lower() for row in rows)),
            'bins':rows}


def upper_envelope(lines, end):
    ordered = sorted(lines,key=lambda row:row[0].mid(),reverse=True)
    assert all(x[0]>y[0] for x,y in zip(ordered,ordered[1:]))
    result = []
    for b,a,i in ordered:
        while result:
            p,q,_,start = result[-1]
            w = (q-a)/(p-b)
            if start is None or w > start:
                break
            assert w < start
            result.pop()
        else:
            w = None
        result.append((b,a,i,w))
    assert all(0<event<end for *_,event in result[1:])
    return result


def integral(profile, end):
    events = [arb(0)]+[w for *_,w in profile[1:]]+[end]
    return sum((r-l)*(a-b*(r+l)/2) for (b,a,_,_),l,r in zip(profile,events[:-1],events[1:]))


def added_area(old,new):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_polyline_obstruction import clip_geometry

    z,D,normal,L = (new[k] for k in ('z','D','normal','length'))
    ar = new['A']
    p2 = ar('0.410')
    prior = coefficients(old)
    all_lines = coefficients(new)
    for b,a,_ in all_lines[len(prior):]:
        # This makes the old envelope unchanged outside the straight segment.
        assert prior[0][0]>b>prior[1][0] and a<prior[0][1] and a<p2*b
    before = upper_envelope(prior,p2)
    after = upper_envelope(all_lines,p2)
    assert before[0][2]==after[0][2]==0 and before[-1][2]==after[-1][2]==1
    lo,hi = after[1][3],after[-1][3]
    knots = [row[3] for row in after[1:]]
    for *_,w in before[1:]:
        assert w<lo or w>hi or lo<w<hi
        if lo<w<hi:
            knots.append(w)

    def point(w,data):
        candidates = [a-b*w for b,a,_ in data]
        y = max(v.lower() for v in candidates).union(max(v.upper() for v in candidates))
        return g.plus(g.plus(z,g.scale(w,D)),g.scale(L*y,normal))

    g.SIGMA = ar('1.5497')*g.PI/180
    inner,_ = polygon('0.0005522',32)
    directions,offsets,*_ = clip_geometry(['0.000555','0.0005522','0.0075'])
    inner.extend((f'R126_EH_{i}',u,h) for i,(u,h) in enumerate(zip(directions,offsets)))
    points = [point(w,lines) for w in knots for lines in (prior,all_lines)]
    clearance = min((h-g.dot(u,x)).lower() for _,u,h in inner for x in points)
    assert clearance > 0
    gain = L**2*(integral(after,p2)-integral(before,p2))
    assert gain > 0
    return {'old_active_lines':[i for _,_,i,_ in before],
            'new_active_lines':[i for _,_,i,_ in after],
            'new_breakpoints_in_w':[str(w) for *_,w in after[1:]],
            'curve_clearance_lower':str(clearance),
            'area_gain_exact':[str(dyadic(v)) for v in (gain.lower(),gain.upper())]},gain,inner


def run(archive, extra):
    sys.path.insert(0,str(archive))
    old, source, instrumented = sweep_r167(archive)
    new, specialized = new_geometry(source,old,extra)
    margins = transfer(old,new)
    ctx.prec = 160
    shape,gain,inner = added_area(old,new)
    prior_path = archive/'ah-three-line-cut-certificate.json'
    prior = json.loads(prior_path.read_bytes())
    assert Fraction(prior['new_AH_third_pivot_exact'])==R167_PIVOT
    assert Fraction(prior['new_AH_third_angle_exact'])==R167_ANGLE
    a,b = map(Fraction,prior['area_exact_interval'])
    low,high = a-dyadic(gain.upper()),b-dyadic(gain.lower())
    assert low < high < a
    return {'additional_AH_lines_exact':[[str(p),str(d)] for p,d in extra],
            'one_sweep_placement':margins,'shape':shape,
            'area_exact_interval':[str(low),str(high)],'area_decimal_interval':format_area(low,high),
            'R150_R156_witnesses_at_identity':fit_archived_witnesses(archive,new,inner),
            'R167_area_sha256':sha256(prior_path.read_bytes()).hexdigest(),
            'instrumented_sweep_sha256':sha256(instrumented.encode()).hexdigest(),
            'derived_geometry_sha256':sha256(specialized.encode()).hexdigest()}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('extra_lines',type=Path)
    parser.add_argument('output',type=Path)
    opt = parser.parse_args()
    raw_lines = opt.extra_lines.read_bytes()
    extra = [tuple(map(Fraction,line)) for line in json.loads(raw_lines)]
    assert extra and all(0<p<Fraction('0.410') and 0<a<Fraction('0.00060') for p,a in extra)
    certificate = run(opt.archive.resolve(),extra)
    certificate['line_input_sha256'] = sha256(raw_lines).hexdigest()
    certificate['script_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    opt.output.write_text(json.dumps(certificate,indent=2)+'\n')
