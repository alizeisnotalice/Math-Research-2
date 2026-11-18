"""Compute full physical Gibbs, AH and eleven-facet EH areas at a single angle; no placement assumption is used."""
import argparse,json,sys,subprocess,io,contextlib
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import ctx,arb
from certify_ah_margin_envelope import replace_line
from certify_ah_three_line_cut import format_area


def cap(archive,sigma,p,depth,end):
    source=archive/'certify_ah_combined_area.py';code=source.read_text().split('baseline = here / ')[0]
    for a,b in (('sigma = "1.5497"',f'sigma = "{sigma}"'),('pivot_parameter = F("0.38")',f'pivot_parameter = F("{p}")'),('angle = F("0.00033")',f'angle = F("{depth}")'),('start, finish = F("0.037"), F("0.038")',f'start, finish = F(0), F("{end}")')):
        code=replace_line(code,a,b)
    env={'__file__':str(source),'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(code,str(source),'exec'),env)
    gain=env['cap_area'];return gain.lo,gain.hi


def run(archive,owner,inp,out):
    sys.path.append(str(owner))
    from certify_eh_upper_chain import front,b
    from certify_eh_joint_parameters import specialize,gain_endpoint
    cfg=json.loads(inp.read_bytes());sigma=cfg['sigma_degrees'];s=F(sigma);ctx.prec=160
    area=json.loads(subprocess.check_output([sys.executable,'-B',str(archive/'certify_gibbs_area_refined.py'),sigma,'16','64']))
    base0,base1=map(F,area['area_interval'])
    old=json.loads((owner/'eh-two-line-cover-certificate.json').read_bytes());params=json.loads((owner/'eh-upper-chain-parameters.json').read_bytes());_,special=specialize(s)
    assert b.g.SIGMA.overlaps(b.arbq(s)*b.g.PI/180)
    prev=b.eh.clip_geometry(old['parameters'])[2];plus=front(old,params)['gain'];eh=prev+plus
    el,eu=map(gain_endpoint,(eh.lower(),eh.upper()))
    rows=[]
    for pivot,depth,end in cfg['AH_chords']:
        lo,hi=cap(archive,sigma,pivot,depth,end);L,U=base0-el-hi,base1-eu-lo
        # Lower endpoint subtracts the two upper gain endpoints.
        L,U=base0-eu-hi,base1-el-lo;assert L<U
        rows.append({'AH_chord':[pivot,depth],'AH_gain_exact_interval':[str(lo),str(hi)],'AH_gain_decimal_interval':format_area(lo,hi),'full_area_exact_interval':[str(L),str(U)],'full_area_decimal_interval':format_area(L,U)})
        print(pivot,depth,rows[-1]['full_area_decimal_interval'],flush=True)
    doc={'parameters':cfg,'gibbs_base':area,'gibbs_base_decimal_interval':format_area(base0,base1),'EH_gain_exact_interval':[str(el),str(eu)],'EH_gain_decimal_interval':format_area(el,eu),'EH_facets':params,'EH_model_sha256':special,'rows':rows,'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':b.g.arb.__module__,'area_only':True}
    out.write_text(json.dumps(doc,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('owner',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.owner.resolve(),opt.input,opt.output)
