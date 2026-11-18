"""Uniform coverage and exact area brackets for an ordered family of two AH cuts."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys
from types import FunctionType

from flint import arb, ctx


def read(path):
    return json.loads(path.read_bytes())


def run(archive, input_path, output):
    from certify_ah_physical_floor import geometry
    from certify_ah_single_contact_chain import validate
    from certify_ah_projected_jets import bracket
    from certify_ah_deep_single_line import green as old_green
    from certify_ah_margin_envelope import replace_line
    from certify_ah_three_line_cut import format_area
    from certify_ah_wider_dual_chord import area_and_witnesses

    cfg=read(input_path)
    (tail,front)=cfg['cuts']
    (p0,p1),(d0,d1)=map(lambda row:tuple(map(F,row)),(tail['pivot'],tail['angle']))
    (q0,q1),(e0,e1)=map(lambda row:tuple(map(F,row)),(front['pivot'],front['angle']))
    assert 0<q0<q1<p0<p1<F('.43')
    assert 0<d0<d1<e0<e1<F('.006')
    far,negative,eps=map(F,cfg['physical_parameters'])
    labelled=cfg.get('support_method')=='first_hit'
    if labelled:
        from certify_ah_first_hit_chain import validate
    strong=[(p1,d1),(q1,e1)]
    weak=[(p0,d0),(q0,e0)]

    geom,source,_=geometry(archive,strong,far=far,negative=negative,eps=eps,labelled=labelled)
    A=geom['A']
    intercept_gap=A(q0)*A(e0).tan()-A(p1)*A(d1).tan()
    assert intercept_gap>0
    slope_ratio=min((A(lo).tan()/A(hi).tan()).lower() for lo,hi in ((d0,d1),(e0,e1)))
    assert 0<slope_ratio<=1

    active={j for j,E in enumerate(geom['deficit_by_bin']) if not E.is_zero()}
    seen=set()
    common=[]
    for low,high,l,r,turn,sgn,orient in cfg['classes']:
        for j in range(low,high):
            assert j in active and j not in seen
            stats=validate(geom,j,F(l),F(r),turn,sgn,orient,angular=12,depth=4)
            assert arb(stats['placement_margin_lower'])>0 and arb(stats['root_margin_lower'])>0
            common.append({'bin':j,**stats})
            seen.add(j)
            print('chain',j,stats['placement_margin_lower'],flush=True)
    remaining,verifier=bracket(geom,active-seen)
    assert {row['bin'] for row in remaining}==active-seen
    bins=sorted(common+remaining,key=lambda row:row['bin'])
    strong_margin=min(arb(row['placement_margin_lower']).lower() for row in bins)

    script=replace_line(inspect.getsource(old_green),'F("0.07")',f'F("{far}")')
    env=dict(old_green.__globals__)
    exec(compile(script,'family_green','exec'),env)
    area=FunctionType(area_and_witnesses.__code__,{**area_and_witnesses.__globals__,'green':env['green']})
    ctx.prec=160
    sl,su,ss,witnesses=area(archive,geom,strong)
    weak_geom,weak_source,_=geometry(archive,weak,far=far,negative=negative,eps=eps,labelled=labelled)
    ctx.prec=160
    wl,wu,ws,_=area(archive,weak_geom,weak)
    if 'reference_single_chord' in cfg:
        r,s=map(F,cfg['reference_single_chord'])
        _,removed_upper,_=env['green'](archive,r,s)
        g_lower=F(read(archive/'gibbs-1p5497-refined-18-128.json')['area_interval'][0])
        eh_upper=F(read(archive/'eh-two-line-cover-certificate.json')['summary']['cut_gain_area_interval_exact'][1])
        baseline=g_lower-eh_upper-removed_upper
    else:
        baseline=F(cfg['reference_AH_area_lower'])
    assert sl<su<wl<wu<baseline

    factor=(strong_margin*slope_ratio).lower()
    assert factor>0
    def digest(path):return sha256(path.read_bytes()).hexdigest()
    result={'complete':True,'lines':[[str(p),str(d)] for p,d in strong],
            'family':cfg,'strong_member_area_exact':list(map(str,(sl,su))),
            'strong_member_area_decimal':format_area(sl,su),
            'weak_member_area_exact':list(map(str,(wl,wu))),'weak_member_area_decimal':format_area(wl,wu),
            'all_members_area_bounds_decimal':format_area(sl,wu),'strict_gain_all_members':str(baseline-wu),
            'area_exact_interval':[str(sl),str(su)],'area_decimal_interval':format_area(sl,su),
            'uniform_cut_line_margin_lower':str(factor),'strong_placement_margin_lower':str(strong_margin),
            'normalized_slope_ratio_lower':str(slope_ratio),'uniform_intercept_gap_lower':str(intercept_gap.lower()),
            'strong_shape':ss,'weak_shape':ws,'bins':bins,
            'R150_R156_witnesses_at_identity':witnesses,'input_sha256':digest(input_path),
            'strong_geometry_sha256':sha256(source.encode()).hexdigest(),
            'weak_geometry_sha256':sha256(weak_source.encode()).hexdigest(),
            'regular_verifier_sha256':sha256(verifier.encode()).hexdigest(),
            'source_hashes':{n:digest(archive/n) for n in
                             ('eh-two-line-cover-certificate.json','certify_ah_centered_deficit.py',
                              'certify_ah_single_contact_chain.py','certify_ah_combined_area.py')},
            'python':sys.version,'flint':geom['flint_version']}
    output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','input','output'):
        cli.add_argument(name,type=Path)
    args=cli.parse_args()
    archive=args.archive.resolve()
    sys.path.insert(0,str(archive))
    sys.path.append(str(Path(__file__).resolve().parent.parent))
    run(archive,args.input,args.output)
