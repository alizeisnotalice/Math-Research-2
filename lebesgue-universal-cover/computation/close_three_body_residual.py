"""Close the finite exact residual list of the three-body seed replay."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from itertools import combinations
import json
from math import floor
from pathlib import Path

from flint import arb
import numpy as np
from certify_three_body_refinement import CHECKER_SHA,SEED_SHA,ROOT_LO,ROOT_HI,corner_bounds,split
from check_three_body_fan_dual import pair_max,SCALE


def close(source,stage):
    file=Path(source)/'checker/lbcheck.py';assert sha256(file.read_bytes()).hexdigest()==CHECKER_SHA
    spec=spec_from_file_location('tomar',file);lib=module_from_spec(spec);spec.loader.exec_module(lib)
    raw=stage.read_bytes();prev=json.loads(raw);target=Q(prev['target']);tau=arb(prev['target'])
    assert prev['checker_sha256']==CHECKER_SHA and prev['seed_sha256']==SEED_SHA and prev['root']==[[str(x) for x in row] for row in [ROOT_LO,ROOT_HI]]
    assert [p['part'] for p in prev['parts']]==list(range(16))
    todo=[(p['part'],i,(tuple(map(Q,a)),tuple(map(Q,b)))) for p in prev['parts'] for i,(a,b) in enumerate(p['unresolved'])]
    assert len(todo)==prev['unresolved'] and all(all(ra<=a<b<=rb for a,b,ra,rb in zip(l,r,ROOT_LO,ROOT_HI)) for _,_,(l,r) in todo)
    checker=lib.Checker([3,5],float(target),256,1e-5);checker.tau=tau
    shifts=[(0,0,0,0,0)]+[tuple(s if i==j else 0 for i in range(5)) for j in range(5) for s in (-1,1)]
    def check(l,r,depth):
        c=[(a+b)/2 for a,b in zip(l,r)];h=[(b-a)/2 for a,b in zip(l,r)]
        assert all(Q(float(x))==x for x in c+h)
        c=np.array(c,dtype=float);h=np.array(h,dtype=float);rows={}
        for shift in shifts:
            ss=checker.verified_runs(c,h,c+h*np.array(shift))
            if ss[0] is None:continue
            row=[]
            for lower in corner_bounds(*ss,h):
                value=floor(float(lower.lower())*SCALE)-1;assert lower>arb(value)/SCALE;row.append(value)
            if min(row)>target*SCALE:
                return {'method':'single','shift':shift,'bound':str(Q(min(row),SCALE))}
            for other,a in rows.items():
                val,w=pair_max(a,row)
                if val>target*SCALE:
                    assert min((w*x+(1-w)*y)/SCALE for x,y in zip(a,row))==val/SCALE
                    return {'method':'dual','shifts':[other,shift],'weight_first':str(w),'bound':str(val/SCALE)}
            rows[shift]=row
        assert depth<8,'unresolved after second stage'
        axis=max(range(5),key=lambda i:h[i]*(0.22 if i==2 else 1))
        return {'method':'split','axis':axis,'children':[check(ll,rr,depth+1) for ll,rr in split(l,r,axis)]}
    results=[{'part':p,'residual_index':i,'proof':check(*box,0)} for p,i,box in todo]
    return {'source_stage_sha256':sha256(raw).hexdigest(),'target':str(target),'residual_boxes':len(todo),
            'completed':True,'proofs':results,'scope':'M(D,B3,B5) >= target, in conjunction with fresh first-stage replay and round-35 reduction'}

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('source');ap.add_argument('stage',type=Path);args=ap.parse_args();print(json.dumps(close(args.source,args.stage),indent=2))
