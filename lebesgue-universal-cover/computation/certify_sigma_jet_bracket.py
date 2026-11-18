"""Test an actual transformed AH root bracket using pinned supports and oriented averaged chart derivatives."""
import argparse,inspect,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import arb
from certify_sigma_extended_graph import geometry
from certify_sigma_pinned_fan import combined,conditional
from certify_ah_local_physical_support import validate as ordinary
from certify_sigma_reachable_roots import reach


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());p,d=map(F,cfg['line']);far,negative,eps=map(F,cfg['physical_parameters']);ns,_,_=geometry(archive,[(p,d)],sigma=cfg['sigma_degrees'],far=far,negative=negative,eps=eps);ns['_sharpen_jets']=True
    env={**ordinary.__globals__,'support_bound':combined};exec(compile(inspect.getsource(ordinary),'pinned_oriented_root','exec'),env);direct=env['validate'];rows=[]
    for j,mode,params in cfg['bins']:
        subsets=reach(ns,j);used=[]
        for a,b in subsets:
            if mode=='direct':
                l,r,t,s,sg=params;used.append(direct(ns,j,F(l),F(r),t,s,sg,begin=a,end=b,angular=12,depth=4))
            else:
                l,anchor,other,r,tau=params
                for low,high,turn,sign in [(l,anchor,300,-1),(other,r,180,1)]:
                    used.append(conditional(ns,j,F(low),F(high),turn,1,1,begin=a,end=b,angular=12,depth=4,anchor=F(anchor),threshold=F(tau),condition_sign=sign))
        row={'bin':j,'reachable':subsets,'classes':used,'positive':all(arb(q[key])>0 for q in used for key in ['root_margin_lower','placement_margin_lower'])};rows.append(row);print(j,row['positive'],[(q['root_margin_lower'],q['placement_margin_lower'],q.get('placement_id')) for q in used],flush=True)
    obj={'parameters':cfg,'rows':rows,'all_passed':all(r['positive'] for r in rows),'source_hashes':{name:sha256((archive/name).read_bytes()).hexdigest() for name in ['certify_sigma_pinned_fan.py','certify_sigma_jet_bracket.py','certify_sigma_extended_graph.py']},'flint':ns['flint_version'],'bits':128,'python':sys.version};out.write_text(json.dumps(obj,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);o=ap.parse_args();run(o.archive.resolve(),o.input,o.output)
