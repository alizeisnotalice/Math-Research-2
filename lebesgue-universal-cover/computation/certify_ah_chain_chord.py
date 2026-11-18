"""Validate a complete single chord with graph-height buffers and the six-variable contact simplex."""
import argparse,contextlib,inspect,io,json,sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb,ctx
from certify_ah_graph_endpoint import geometry as endpoint_geometry
from certify_ah_single_contact_chain import validate
from certify_ah_projected_jets import bracket
from certify_ah_margin_envelope import coefficients,replace_line
from certify_ah_three_line_cut import fit_archived_witnesses,format_area
import certify_ah_deep_single_line as base


def geometry(archive,extra,*,far,negative,eps):
    template=inspect.getsource(endpoint_geometry)
    marker='    original=replace_line(original,marker,extra_code+marker);'
    body="    extra_code+='    prefix=replace_line(prefix,'+repr('((-eps,0),(0,eps))')+','+repr('tuple((eps*F(i,2),eps*F(i+1,2)) for i in range(-2,2))')+')\\n'\n"
    # Both coordinate loops use the same two-cell source string, so replace it twice after checking the cardinality.
    body="    extra_code+=\"    assert prefix.count('((-eps,0),(0,eps))')==2\\n    prefix=prefix.replace('((-eps,0),(0,eps))','tuple((eps*F(i,2),eps*F(i+1,2)) for i in range(-2,2))')\\n    prefix=replace_line(prefix,'A(\\\"0.000054\\\")','A(\\\"0.000080\\\")')\\n\"\n"
    template=replace_line(template,marker,body+marker);env=dict(endpoint_geometry.__globals__);exec(compile(template,'subdivided_derivative_chart','exec'),env);return env['geometry'](archive,extra,far=far,negative=negative,eps=eps)


def run(archive,inp,out):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_polyline_obstruction import clip_geometry
    cfg=json.loads(inp.read_bytes());(p,d),=[tuple(map(F,l)) for l in cfg['lines']];far,negative,eps=map(F,cfg['physical_parameters']);ns,full,_=geometry(archive,[(p,d)],far=far,negative=negative,eps=eps)
    active={j for j,E in enumerate(ns['deficit_by_bin']) if not E.is_zero()};common=[];seen=set()
    for low,high,l,r,turn,sgn,orient in cfg['classes']:
        for j in range(low,high):
            assert j in active and j not in seen;seen.add(j);row=validate(ns,j,F(l),F(r),turn,sgn,orient,angular=12,depth=4);assert arb(row['placement_margin_lower'])>0 and arb(row['root_margin_lower'])>0,row;common.append(row);print('joint',j,row['placed_leaves'],row['root_free_leaves'],flush=True)
    regular,regular_src=bracket(ns,sorted(active-seen));assert seen|{x['bin'] for x in regular}==active
    ctx.prec=160;cs=coefficients(ns);s,a,_=cs[-1];A=ns['A'];assert all(a-s*w>h-v*w for v,h,_ in cs[:-1] for w in (-A(negative),A(p)))
    g.SIGMA=A('1.5497')*g.PI/180;inner,_=polygon('0.0005522',32);directions,offsets,*_=clip_geometry(['0.000555','0.0005522','0.0075']);inner.extend((f'R126_EH_{i}',u,h) for i,(u,h) in enumerate(zip(directions,offsets)));guards=fit_archived_witnesses(archive,ns,inner)
    code=inspect.getsource(base.green);code=replace_line(code,'F("0.07")',f'F("{far}")');gn=dict(base.green.__globals__);exec(compile(code,'paired_height_green_area','exec'),gn);L,U,shape=gn['green'](archive,p,d)
    gib=json.loads((archive/'gibbs-1p5497-refined-18-128.json').read_bytes());eh=json.loads((archive/'eh-two-line-cover-certificate.json').read_bytes())['summary'];assert F(eh['sigma_degrees'])==F(gib['sigma_degrees'])==F('1.5497') and F(eh['delta'])==F('.000555') and eh['AH_relocation_after_root'] and eh['delta_is_enclosing_placement_depth']; gl,gu=map(F,gib['area_interval']);el,eu=map(F,eh['cut_gain_area_interval_exact']);lo,hi=gl-eu-U,gu-el-L;assert lo<hi<F('.844093088349096210409615928214643078085')
    data={'parameters':cfg,'physical_parameters':cfg['physical_parameters'],'lines':cfg['lines'],'common_classes':common,'bins':regular,'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),'shape':shape,'R150_R156_witnesses_at_identity':guards,'input_sha256':sha256(inp.read_bytes()).hexdigest(),'geometry_sha256':sha256(full.encode()).hexdigest(),'regular_scan_sha256':sha256(regular_src.encode()).hexdigest(),'source_hashes':{q:sha256((archive/q).read_bytes()).hexdigest() for q in ['certify_ah_single_contact_chain.py','certify_ah_attainable_root.py','certify_ah_paired_height.py','certify_ah_graph_endpoint.py','certify_gibbs_area_refined.py','gibbs-1p5497-refined-18-128.json','eh-two-line-cover-certificate.json']},'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':ns['flint_version'],'bits':128};out.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
