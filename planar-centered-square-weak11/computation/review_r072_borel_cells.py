#!/usr/bin/env python3
"""Independent Borel saved-cell review: literal 2D atoms and exact partition."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ov(x,r,lo,hi):return max(0,min(x+r,hi)-max(x-r,lo))


def main():
    ap=argparse.ArgumentParser()
    for k in ('summary','certificate','producer','output'):ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--mode',choices=['all','sample'],required=True);ap.add_argument('--sample-per-state',type=int,default=16);a=ap.parse_args();start=time.perf_counter()
    s=json.loads(a.summary.read_text());raw=json.loads(gzip.decompress(a.certificate.read_bytes()));assert sha(a.certificate)==s['certificate_sha256'] and sha(a.producer)==s['receipt']['script_sha256']
    n=s['n'];D=s['integer_coordinate_scale'];domain=tuple(raw['receiver_domain_units']);J=tuple(raw['J_axis_units']);grid=raw['grid_units'];points=list(product(grid,grid));assert len(points)==n*n and raw['n']==n
    assert grid==[i*(n+1)*D//(n*n) for i in range(n)] and J==(-D//(n*n),D)
    assert domain==(-D//2,(n*n-1)*D//(n*n)+D//2)*2
    trie={};areas=Counter();counts=Counter();groups={-1:[],0:[],1:[]}
    for row in raw['leaves']:
        lo,hi,bot,top=row['cell_units'];depth=row['depth'];state=row['state'];assert lo<hi and bot<top and state in groups
        div=1<<depth;wx=(domain[1]-domain[0])//div;wy=(domain[3]-domain[2])//div
        assert wx*div==domain[1]-domain[0] and wy*div==domain[3]-domain[2]
        ix=(lo-domain[0])//wx;iy=(bot-domain[2])//wy
        assert 0<=ix<div and 0<=iy<div and (lo,hi)==(domain[0]+ix*wx,domain[0]+(ix+1)*wx) and (bot,top)==(domain[2]+iy*wy,domain[2]+(iy+1)*wy)
        node=trie
        for bit in range(depth-1,-1,-1):
            assert 'leaf' not in node
            child=((ix>>bit)&1)+2*((iy>>bit)&1);node=node.setdefault(child,{})
        assert not node;node['leaf']=True
        areas[state]+=(hi-lo)*(top-bot);counts[state]+=1;groups[state].append(row)
    def complete(node):
        if 'leaf' in node:assert len(node)==1
        else:
            assert set(node)=={0,1,2,3}
            for child in node.values():complete(child)
    complete(trie)
    assert sum(areas.values())==(domain[1]-domain[0])*(domain[3]-domain[2])
    assert {str(k):counts[k] for k in (-1,0,1)}==s['leaf_counts']
    assert F(areas[1],D*D)==F(s['area_lower']) and F(areas[1]+areas[0],D*D)==F(s['area_upper']) and F(areas[0],D*D)==F(s['unresolved_area'])
    checked={-1:0,1:0};literal_atom_distance_tests=0;outer_event_tests=0;mid=(J[0]+J[1])//2;assert 2*mid==sum(J)
    for state in (-1,1):
        records=groups[state]
        if a.mode=='sample' and len(records)>a.sample_per_state:
            ids=sorted({i*(len(records)-1)//(a.sample_per_state-1) for i in range(a.sample_per_state)});records=[records[i] for i in ids]
        for row in records:
            u,v,b,t=row['cell_units'];near=(min(v,max(u,mid)),min(t,max(b,mid)));far=(u if abs(u-mid)>=abs(v-mid) else v,b if abs(b-mid)>=abs(t-mid) else t)
            if state==1:
                R=row['witness']['radius_units'];assert 0<=2*R<D
                count=sum(max(abs(px-u),abs(px-v),abs(py-b),abs(py-t))<=R for px,py in points);literal_atom_distance_tests+=len(points)
                overlap=ov(near[0],R,*J)*ov(near[1],R,*J)
                margin=count*(n+1)*D*D-n**3*overlap-4*n*n*R*R
                assert margin>0 and count==row['witness']['guaranteed_atom_count'] and overlap==row['witness']['maximum_J_overlap_units_squared'] and margin==row['witness']['strict_integer_margin']
            else:
                # Actual 2D source-to-cell minimum distances; this avoids the
                # producer's two separately sorted 1D possible-distance lists.
                distances=Counter(max(u-px,px-v,b-py,py-t,0) for px,py in points);literal_atom_distance_tests+=len(points)
                count=0;minimum=None
                for R,multiplicity in sorted(distances.items()):
                    if 2*R>=D:break
                    count+=multiplicity;overlap=ov(far[0],R,*J)*ov(far[1],R,*J)
                    margin=n**3*overlap+4*n*n*R*R-count*(n+1)*D*D;assert margin>=0
                    minimum=margin if minimum is None else min(minimum,margin);outer_event_tests+=1
                assert minimum==row['witness']['minimum_nonnegative_integer_margin']
            checked[state]+=1
    fully=a.mode=='all'
    result=dict(status='PASS_INDEPENDENT_BOREL_COMPLETE_PARTITION_AND_LEAF_PREDICATES' if fully else 'PASS_INDEPENDENT_BOREL_COMPLETE_PARTITION_AND_SAMPLED_LEAF_PREDICATES',source_sha256=sha(Path(__file__)),python=sys.version,summary_sha256=sha(a.summary),certificate_sha256=sha(a.certificate),producer_sha256=sha(a.producer),n=n,all_partition_leaves=len(raw['leaves']),all_leaf_counts={str(k):counts[k] for k in (-1,0,1)},independently_replayed_classified_leaves={str(k):v for k,v in checked.items()},literal_2D_atom_distance_tests=literal_atom_distance_tests,outer_actual_2D_distance_event_tests=outer_event_tests,certified_area_lower=s['area_lower'],certified_area_upper=s['area_upper'],independently_full_area_bracket_verified=fully,n5_atomic_lim_upper='1927/1250' if n==5 else None,n5_strict_discontinuity_gap_lower=str(F(s['area_lower'])-F(1927,1250)) if n==5 and fully else None,elapsed_seconds=time.perf_counter()-start,scope='All leaf partitions and exact area bookkeeping checked. '+('Every classified leaf independently replayed by literal 2D atom counts; full bracket independently verified.' if fully else 'Only deterministic quantile sample of each classified state replayed; no independent full bracket claim, unreplayed classification relies on producer proof/code/results.')+' Uses separately reviewed global Borel cap/Hresidual-empty lemma; not a finite-atomic outside counterexample or a general C2 bound.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],n=n,partition_leaves=result['all_partition_leaves'],classifications_checked=checked,seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
