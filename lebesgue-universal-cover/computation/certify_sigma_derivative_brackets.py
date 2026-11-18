"""Verify the high/low root split selected by the real initial compatibility derivative."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import arb
from certify_sigma_extended_graph import geometry
from certify_sigma_derivative_classes import validate
from certify_sigma_reachable_roots import reach


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());p,d=map(F,cfg['line']);far,negative,eps=map(F,cfg['physical_parameters']);n,src,_=geometry(archive,[(p,d)],sigma=cfg['sigma_degrees'],far=far,negative=negative,eps=eps);n['_sharpen_jets']=True;rows=[]
    for j,pr in cfg['slope_bins']:
        cut,l,r,L,R=map(F,pr);used=[];ranges=reach(n,j)
        for x,y in ranges:
            for lo,hi,t,s,kind,side in [(l,r,0,1,-1,-1),(L,R,120,-1,1,1)]:
                q=validate(n,j,lo,hi,t,s,kind,cut=cut,side=side,begin=x,end=y,angular=12,depth=4);q['initial_slope_class']={'comparison':side,'cut':str(cut)};q.pop('anchor_degrees');q.pop('anchor_comparison_sign');q.pop('anchor_threshold');used.append(q)
        row={'bin':j,'reachable':ranges,'classes':used,'positive':all(arb(q[k])>0 for q in used for k in ['root_margin_lower','placement_margin_lower'])};rows.append(row);print(j,row['positive'],[(q['root_margin_lower'],q['placement_margin_lower'],q['placement_id']) for q in used],flush=True)
    result={'parameters':cfg,'rows':rows,'all_passed':all(r['positive'] for r in rows),'geometry_sha256':sha256(src.encode()).hexdigest(),'source_hashes':{k:sha256((archive/k).read_bytes()).hexdigest() for k in ['certify_sigma_derivative_classes.py','certify_sigma_derivative_brackets.py','certify_sigma_pinned_fan.py','certify_sigma_extended_graph.py']},'python':sys.version,'flint':n['flint_version'],'bits':128};out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('archive',type=Path);a.add_argument('input',type=Path);a.add_argument('output',type=Path);o=a.parse_args();run(o.archive.resolve(),o.input,o.output)
