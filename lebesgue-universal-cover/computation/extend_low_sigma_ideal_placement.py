"""Extend the low Gibbs angle ideal-body theorem past the proposed AH pivot."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import certify_two_sigma_ideal_placement as base

PARTS = [
    ('.42849','.42875','53.46','53.56','54.13','54.22'),
    ('.42875','.429','53.51','53.62','54.08','54.18'),
    ('.429','.42925','53.57','53.70','54.01','54.13'),
    ('.42925','.4294','53.64','53.76','53.95','54.06'),
    ('.4294','.42949','53.71','53.815','53.892','54'),
    ('.42949','.429515','53.77','53.842','53.865','53.94'),
]


def certify():
    sigma=F('1.54')
    data=base.xy(sigma)
    D=data[5]
    rows,_,_,_=base.edges(data)
    results=[]
    last=F(PARTS[0][0])
    for start,stop,x1,x2,y1,y2 in PARTS:
        a,b=map(F,(start,stop))
        assert a==last and a<b and F(x2)<F(y1)
        last=b
        branch=[]
        for k,(low,high) in enumerate(((x1,x2),(y1,y2))):
            l,r=map(F,(low,high))
            probes=(l,r+60) if k==0 else (l+60,r)
            root_cells=128
            sign=[]
            for i in range(root_cells):
                lo,hi=(a+(b-a)*F(j,root_cells) for j in (i,i+1))
                sign.extend(base.rig(sigma,sigma,lo,hi,probe_angles=probes))
            assert min(sign)>0,(start,stop,low,high,str(min(sign)))
            placement_cells=32
            margins=[]
            for i in range(placement_cells):
                lo,hi=(a+(b-a)*F(j,placement_cells) for j in (i,i+1))
                margins.append(base.validate_cell(lo,hi,l,r,D,rows))
            lower,name=min(margins)
            assert lower>0,(start,stop,low,high,str(lower),name)
            branch.append({'angles':[low,high],'sign_lower':str(min(sign)),
                           'minimum_slack_lower':str(lower),'wall':name,
                           'root_cells':root_cells,'placement_cells':placement_cells})
        results.append({'w_interval':[start,stop],'branches':branch})
    return results

if __name__=='__main__':
    output=Path(sys.argv[1])
    output.write_text(json.dumps({'complete':True,'sigma':'1.54','AH_pivot':'.429','AH_depth':'.009',
                                  'EH_facets':11,'cells':certify(),
                                  'base_script_sha256':hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
                                  'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
