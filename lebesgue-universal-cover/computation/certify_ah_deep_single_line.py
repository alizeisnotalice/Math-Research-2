"""Certify one AH chord across the paired arc and straight segment, with directed error and true Green area."""
import argparse
import contextlib
from fractions import Fraction as F
from hashlib import sha256
import io,json,sys
from pathlib import Path
from flint import ctx
from certify_ah_fine_pivot import make_geometry, fine_verify
from certify_ah_margin_envelope import coefficients,replace_line
from certify_ah_three_line_cut import fit_archived_witnesses,format_area


def green(path,p,angle):
    source=path/'certify_ah_combined_area.py'; code=source.read_text().split('baseline = here / ')[0]
    code=replace_line(code,'pivot_parameter = F("0.38")',f'pivot_parameter = F("{p}")'); code=replace_line(code,'angle = F("0.00033")',f'angle = F("{angle}")'); code=replace_line(code,'start, finish = F("0.037"), F("0.038")','start, finish = F(0), F("0.07")'); ns={'__file__':str(source),'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(code,str(source),'exec'),ns)
    cap=ns['cap_area'];return cap.lo,cap.hi,{'crossing_parameter_exact':[str(ns['start']),str(ns['finish'])],'crossing_derivative_lower':str(ns['slope'].lo),'Taylor_order':ns['order'],'Taylor_cells':ns['cells'],'arc_remainder_absolute_upper':str(ns['remainder_abs']),'cap_area_exact_interval':[str(cap.lo),str(cap.hi)],'green_specialization_sha256':sha256(code.encode()).hexdigest()}


def run(archive,inp,out):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_polyline_obstruction import clip_geometry
    extra=[tuple(map(F,r)) for r in json.loads(inp.read_bytes())];(p,d),=extra; assert 0<d<=F('.0006')
    ns,geo=make_geometry(archive,extra); checked,angles=fine_verify(ns,None)
    ctx.prec=160; cs=coefficients(ns); b,a,_=cs[-1]; A=ns['A']; minw=-A(ns['pbound']); maxw=A(p)
    assert all(a-b*w>old_a-old_b*w for old_b,old_a,_ in cs[:-1] for w in (minw,maxw))
    g.SIGMA=A('1.5497')*g.PI/180; inner,_=polygon('0.0005522',32); dirs,offsets,*_=clip_geometry(['0.000555','0.0005522','0.0075']); inner.extend((f'R126_EH_{i}',u,h) for i,(u,h) in enumerate(zip(dirs,offsets)))
    guards=fit_archived_witnesses(archive,ns,inner); cap_lo,cap_hi,shape=green(archive,p,d); shape['old_AH_lines_dominated_on_physical_range']=[str(-ns['pbound']),str(p)]
    def read(name): return json.loads((archive/name).read_bytes())
    gibbs=read('gibbs-1p5497-refined-18-128.json'); eh=read('eh-two-line-cover-certificate.json')['summary'];assert F(eh['sigma_degrees'])==F(gibbs['sigma_degrees'])==F('1.5497') and F(eh['delta'])==F('.000555') and eh['AH_relocation_after_root'] and eh['delta_is_enclosing_placement_depth']; a,b=map(F,gibbs['area_interval']); el,ehigh=map(F,eh['cut_gain_area_interval_exact']); lo,hi=a-ehigh-cap_hi,b-el-cap_lo; assert lo<hi<F('.844093390960250720226452027837524452015')
    deps=['gibbs-1p5497-refined-18-128.json','eh-two-line-cover-certificate.json','certify_gibbs_area_refined.py','certify_ah_combined_area.py','certify_ah_fine_pivot.py','certify_ah_adaptive_bracket.py']; obj={'lines':json.loads(inp.read_bytes()),'geometry_sha256':sha256(geo.encode()).hexdigest(),'fine_verifier_sha256':sha256(angles.encode()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'source_hashes':{name:sha256((archive/name).read_bytes()).hexdigest() for name in deps},'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':ns['flint_version'],'geometry_and_placement_bits':128,'bins':checked,'shape':shape,'area_exact_interval':[str(lo),str(hi)],'area_decimal_interval':format_area(lo,hi),'R150_R156_witnesses_at_identity':guards};out.write_text(json.dumps(obj,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
