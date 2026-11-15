#!/usr/bin/env python3
"""Root-run exact equivalence gate: optimized geometry/cuts vs frozen R61."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def independent_difference(left,right):
    """Independent endpoint-cell membership enumeration, for small tests only."""
    rs=left+right
    xs=sorted({v for r in rs for v in r[:2]})
    ys=sorted({v for r in rs for v in r[2:4]})
    area=0
    for a,b in zip(xs,xs[1:]):
        for c,d in zip(ys,ys[1:]):
            def inside(rects):
                return any(2*r[0]<a+b<2*r[1] and 2*r[2]<c+d<2*r[3] for r in rects)
            if inside(left) and not inside(right):area+=(b-a)*(d-c)
    return area


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--saved-inputs',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    newpath=Path(__file__).with_name('exact_source_frontier.py')
    new=load('frontier_new',newpath);old=load('frontier_frozen_r061',args.reference)
    saved=json.loads(args.saved_inputs.read_text());cases=[]
    # One non-control saved input in each of the actual twelve independent batches.
    for i,batch in enumerate(saved['batches']):
        cases.append((f'saved-batch-{i}-case-{1+i%7}',batch['inputs'][1+i%7]))
    minimal=json.loads(Path(__file__).with_name('minimal-input.json').read_text())
    cases.append(('new-minimal',minimal['strata'][0]['inputs'][0]))
    start=time.perf_counter();checked=[]
    for name,inp in cases:
        p=inp['points_quarters'];m=inp['mass_sixteenths'];bits=24
        a,b,c=old.rectangles(p,m,bits);x,y,z,gamma=new.rectangles(p,m,bits)
        assert sorted(a)==sorted(x) and sorted(b)==sorted(y) and c==z,name
        wa,sa=old.cells(a,b,len(p));wn,sn=new.cells(x,y,len(p))
        assert wa==wn and sa==sn,name
        ra,za,ma=old.cuts(wa,p,m,bits);rn,zn,mn=new.cuts(wn,p,m,bits)
        assert za==zn and ma==list(mn),name
        assert all(rn[k]==v for k,v in ra.items()),(name,ra,rn)
        # The exact flow is generated and self-checked by the new program.
        flow=new.exact_flow(wn,m,zn,mn,rn['raw_critical_mask'])
        assert F(rn['geometric_bound'])<=F(rn['raw_bound'])
        checked.append({'name':name,'sources':len(p),'closed_masks':c,
                        'all_cuts':rn['all_cuts'],'bound':rn['geometric_bound'],
                        'flow_arcs':len(flow['arcs'])})
    # Tiny direct membership check exercises overlapping/empty difference regions.
    left=[(0,7,0,5,1),(3,9,2,8,2)]
    right=[(1,4,-1,6,4),(6,10,4,9,8)]
    direct=independent_difference(left,right)
    assert new.difference_area(left,right)==direct
    assert new.difference_area(left,left)==0
    payload={'status':'PASS','scope':'Exact equivalence on 12 saved R61 inputs and one new minimal input; independent rectangle-difference check.',
             'cases':checked,'difference_area':direct,'wall_seconds':time.perf_counter()-start,
             'python':sys.version,'files':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in (Path(__file__),newpath,args.reference,args.saved_inputs)}}
    new.atomic_json(args.out,payload);print(json.dumps({'status':'PASS','cases':len(checked),'seconds':payload['wall_seconds']}))


if __name__=='__main__':main()
