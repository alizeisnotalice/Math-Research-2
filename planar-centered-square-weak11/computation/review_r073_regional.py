#!/usr/bin/env python3
"""New nine-region ledger from exact atoms; independent horizontal area merge."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from math import lcm
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    for k in ('input','saved','producer','helper','output'):ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--case',default='contracted-unique');a=ap.parse_args();start=time.perf_counter()
    spec=importlib.util.spec_from_file_location('frozen_independent_r072_geometry',a.helper);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    raw=json.loads(a.input.read_text());packet=json.loads(a.saved.read_text());assert packet['source_sha256']==sha(a.producer);saved=next(c for c in packet['cases'] if c['name']==a.case);assert sha(a.input)==saved['input_sha256']
    P=[tuple(map(F,p)) for p in raw['positions']];m=list(map(F,raw['masses']));assert all(w>0 for w in m) and len(set(P))==len(P)
    S=set(raw['S_indices']);assert len(S)==25 and sum((m[i] for i in S),F(0))==1
    base=lcm(*(c.denominator for p in P for c in p),*(v.denominator for v in m));D=base*(1<<saved['source_geometry_bits']);assert D==saved['coordinate_scale']
    full,scans=helper.direct_closures(P,m);R=[i for i in range(len(P)) if i not in S];residual,rs=helper.direct_closures([P[i] for i in R],[m[i] for i in R]);fi,fo,_=helper.receivers(full,D);ri,ro,_=helper.receivers(residual,D)
    assert max(row['mass'] for row in full if row['eligible'])==1
    endpoints=[F(-1,2),F(-1,25),F(1),F(73,50)];regions=[];comparisons=0
    def clip(rr,b):return [(max(l,b[0]),min(r,b[1]),max(d,b[2]),min(t,b[3])) for l,r,d,t in rr if max(l,b[0])<min(r,b[1]) and max(d,b[2])<min(t,b[3])]
    for ix,iy in product(range(3),repeat=2):
        bb=(endpoints[ix],endpoints[ix+1],endpoints[iy],endpoints[iy+1]);box=tuple(int(x*D) for x in bb);assert all(F(x,D)==y for x,y in zip(box,bb))
        area=(bb[1]-bb[0])*(bb[3]-bb[2]);(flo,rhi,dlo),_=helper.y_slabs(clip(fi,box),clip(ro,box));(fhi,rlo,dhi),_=helper.y_slabs(clip(fo,box),clip(ri,box))
        vals={'full_halo':(F(flo,D*D),F(fhi,D*D)),'residual_halo':(F(rlo,D*D),F(rhi,D*D)),'deletion':(F(dlo,D*D),F(dhi,D*D)),'missing_deletion':(area-F(dhi,D*D),area-F(dlo,D*D)),'full_halo_vacancy':(area-F(fhi,D*D),area-F(flo,D*D))}
        tag=f'x{ix}y{iy}';entry=saved['regions'][tag];assert tuple(map(F,entry['bounds']))==bb and F(entry['area'])==area
        for k,(lo,hi) in vals.items():assert (F(entry[k]['lower']),F(entry[k]['upper']))==(lo,hi);comparisons+=2
        regions.append(dict(id=tag,area=str(area),kind=entry['kind'],metrics={k:dict(lower=str(lo),upper=str(hi)) for k,(lo,hi) in vals.items()}))
    # Every nonempty closed-unit-window capture pattern: an independent
    # arrangement of each source's origin-rectangle projections, using masks.
    patterns=[]
    for axis in (0,1):
        boundaries=sorted({p[axis]-1 for p in P}|{p[axis] for p in P});reps=boundaries+[(x+y)/2 for x,y in zip(boundaries,boundaries[1:])]
        patterns.append({sum(1<<i for i,p in enumerate(P) if p[axis]-1<=z<=p[axis]) for z in reps})
    captures={x&y for x,y in product(*patterns)}-{0};saturating=[];minimum=F(1)
    for mask in captures:
        mass=sum((v for i,v in enumerate(m) if mask>>i&1),F(0));assert mass<=1;minimum=min(minimum,1-mass)
        if mass==1:saturating.append(mask)
    oldwindows=saved['unit_window_caps'];assert len(captures)==oldwindows['distinct_nonempty_capture_masks'];assert set(saturating)=={int(r['capture_mask']) for r in oldwindows['saturating_unit_window_masks']};assert minimum==F(oldwindows['global_minimum_unit_window_slack'])
    result=dict(status='PASS_NEW_ATOM_BASED_NINE_REGION_AND_UNIT_WINDOW_AUDIT',source_sha256=sha(Path(__file__)),python=sys.version,input_sha256=sha(a.input),saved_result_sha256=sha(a.saved),producer_sha256=sha(a.producer),frozen_independent_helper_sha256=sha(a.helper),case=a.case,positive_sources=len(P),direct_full_coordinate_rectangles=scans,full_closures=len(full),actual_residual_closures=len(residual),new_region_endpoint_comparisons=comparisons,regions=regions,all_nonempty_unit_capture_patterns=len(captures),all_saturating_unit_capture_patterns=len(saturating),minimum_closed_unit_window_slack=str(minimum),elapsed_seconds=time.perf_counter()-start,scope='One selected unique38 control only, newly clipped nine-region metrics from raw atoms and independently reconstructed closures/roots, with horizontal yslab merge; no unclippedR72wholearea rerun. Complete finite unit-origin capture patterns including ties/point strata; other three controls source/code bound only. No general outer or corner inequality.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],case=a.case,region_endpoints=comparisons,unit_patterns=len(captures),saturating=len(saturating),seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
