#!/usr/bin/env python3
"""Three cut leaders + three deletion leaders. No production imports/rerun."""
import argparse
from array import array
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def merged(intervals):
    result=[]
    for a,b in sorted(intervals):
        if result and a<=result[-1][1]:result[-1]=(result[-1][0],max(result[-1][1],b))
        else:result.append((a,b))
    return result


def difference(left,right):
    """Rebuilt horizontal strips and interval merges, no global source grid."""
    ys=sorted({r[j] for r in left+right for j in (2,3)})
    result=0
    for a,b in zip(ys,ys[1:]):
        l=merged([(r[0],r[1]) for r in left if r[2]<=a and b<=r[3]])
        r=merged([(r[0],r[1]) for r in right if r[2]<=a and b<=r[3]])
        length=sum(d-c for c,d in l);i=j=0
        while i<len(l) and j<len(r):
            length-=max(0,min(l[i][1],r[j][1])-max(l[i][0],r[j][0]))
            if l[i][1]<r[j][1]:i+=1
            else:j+=1
        assert length>=0;result+=(b-a)*length
    return result


def all_cuts_transposed(weights,p,m,bits,saved):
    """Push each sparse label to all supersets; NOT a zeta transform."""
    n=len(m);N=1<<n;full=N-1;z=[0]*N
    for mask,w in weights.items():
        complement=full^mask;sub=complement
        while True:
            z[mask|sub]+=w
            if not sub:break
            sub=(sub-1)&complement
    mass=array('q',[0])*N
    best_num=0;best_den=1;best_mask=0
    for s in range(1,N):
        bit=s&-s;mass[s]=mass[s^bit]+m[bit.bit_length()-1]
        if z[s]*best_den>best_num*mass[s]:best_num,best_den,best_mask=z[s],mass[s],s
    raw=F(16*best_num,(1<<(2*bits))*best_den)
    assert raw==F(saved['raw_bound']) and best_mask==saved['raw_critical_mask']
    # For these actual leaders raw<B. Every R53 terminal ratio is >=B,
    # so every min equals its raw ratio and no terminal improves any cut.
    assert raw<F(saved.get('B','1994727599733432151936139/604462909807314587353088'))
    assert F(saved['geometric_bound'])==raw and saved['terminal_improved_cuts']==0
    return {'raw':str(raw),'raw_mask':best_mask,'z_at_raw':best_num,
            'm_at_raw':best_den,'all_cuts':N-1}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--reviewer',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    spec=importlib.util.spec_from_file_location('independent_r061',args.reviewer)
    ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
    started=time.perf_counter();cutrows=[];delrows=[];pins=[]
    for name,target in [('allcuts-batch',cutrows),('deletion-batch',delrows)]:
        for path in sorted((args.root/name).glob('stratum-*.json')):
            data=json.loads(path.read_text());assert data['status']=='PASS'
            pins.append({'path':str(path),'sha256':sha(path)})
            target.extend(data['cases'])
    cutrows.sort(key=lambda r:F(r['geometric_bound']),reverse=True)
    delrows=[r for r in delrows if r['deletions']['max_mass_min_upper'] is not None]
    delrows.sort(key=lambda r:F(r['deletions']['max_mass_min_upper']),reverse=True)
    result={'status':'RUNNING','scope':'Three cut leaders and three deletion leaders only; not all72 independent geometry replay.',
            'source_sha256':sha(Path(__file__)),'reviewer_sha256':sha(args.reviewer),'source_output_pins':pins,
            'python':sys.version,'cut_leaders':[],'deletion_leaders':[]}
    for saved in cutrows[:3]:
        p=saved['input']['points_quarters'];m=saved['input']['mass_sixteenths'];bits=saved['bits']
        rects,closed,gamma=ref.boxes(p,m,bits);weights,pairs,ncells=ref.overlay(rects)
        assert closed==saved['closed_source_masks'] and gamma==int(saved['Gamma_numerator'])
        assert weights=={s:int(w) for s,w in saved['area_weight_numerators']}
        cut=all_cuts_transposed(weights,p,m,bits,saved)
        arcs=ref.flow_audit(saved['flow'],weights,m,cut)
        result['cut_leaders'].append({'id':saved['input']['id'],'all_cuts':cut['all_cuts'],
                                      'bound':cut['raw'],'flow_arcs':arcs,'independent_cells':ncells})
    for saved in delrows[:3]:
        p=saved['input']['points_quarters'];m=saved['input']['mass_sixteenths'];bits=saved['bits']
        rects,closed,gamma=ref.boxes(p,m,bits)
        inn=[r[:5] for r in rects if r[5]==0];out=[r[:5] for r in rects if r[5]==1]
        maskmass={r[4]:sum(m[i] for i in range(len(m)) if r[4]>>i&1) for r in out}
        largest=max(maskmass.values());ties=sorted(s for s,v in maskmass.items() if v==largest)
        d=saved['deletions'];assert largest==d['maximum_strictly_eligible_mass_sixteenths']
        assert ties==d['all_maximum_strictly_eligible_masks']
        assert difference(inn,[])==int(saved['full_inner_area_numerator'])
        assert difference(out,[])==int(saved['full_outer_area_numerator'])
        for block in d['blocks']:
            keep=[i for i in range(len(m)) if not block['mask']>>i&1]
            rr,cc,gg=ref.boxes([p[i] for i in keep],[m[i] for i in keep],bits)
            ri=[r[:5] for r in rr if r[5]==0];ro=[r[:5] for r in rr if r[5]==1]
            assert cc==block['remaining_closed_source_masks'] and gg==int(block['remaining_Gamma_numerator'])
            assert difference(ri,[])==int(block['remaining_inner_area_numerator'])
            assert difference(ro,[])==int(block['remaining_outer_area_numerator'])
            assert difference(inn,ro)==int(block['loss_numerator_lower'])
            assert difference(out,ri)==int(block['loss_numerator_upper'])
        result['deletion_leaders'].append({'id':saved['input']['id'],'maximum_masks':ties,
                                           'all_saved_blocks_recomputed':len(d['blocks']),
                                           'minimum_maxmass_upper':d['max_mass_min_upper']})
    result['status']='PASS';result['wall_seconds']=time.perf_counter()-started
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','seconds':result['wall_seconds'],'cut_leaders':len(result['cut_leaders']),'deletion_leaders':len(result['deletion_leaders'])}))


if __name__=='__main__':main()
