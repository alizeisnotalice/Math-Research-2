#!/usr/bin/env python3
"""Separate R66 high-competitor erosion-ring diagnostic, no budget-source edit.

For T: a=m(T intersect S), b=m(T minus S). The full-T bbox receiver
rectangle at mass b is contained in the remaining halo. Thus true unpaid
E_T is contained in the rectangle ring E_T(m) minus E_T(b). Aggregate rings
are computed as a UNION of rectangle differences, never as a sum of a's.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

BUDGET_SHA='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
GEOMETRY_SHA='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def radii(m,bits):
    if m==0:return 0,0
    D=1<<bits;v=m.numerator*D*D;d=4*m.denominator;k=math.isqrt(v//d)
    assert d*k*k<=v<d*(k+1)*(k+1)
    return max(0,k-1),k+1


def receiver(bounds,m,bits,mask):
    a,b,c,d=bounds
    if m==0 or max(b-a,d-c)**2>=m:return [],[]
    D=1<<bits;q=[v*D for v in bounds];assert all(v.denominator==1 for v in q)
    a,b,c,d=map(int,q);lo,hi=radii(m,bits);result=[]
    for r in (lo,hi):
        rect=(b-r,a+r,d-r,c+r,mask)
        result.append([rect] if rect[0]<rect[1] and rect[2]<rect[3] else [])
    return tuple(result)


def subtract_rect(A,B):
    """Disjoint rectangular decomposition of A minus B, up to boundaries."""
    if not B:return [A]
    B=B[0];x0=max(A[0],B[0]);x1=min(A[1],B[1]);y0=max(A[2],B[2]);y1=min(A[3],B[3])
    if x0>=x1 or y0>=y1:return [A]
    candidates=[(A[0],x0,A[2],A[3],A[4]),(x1,A[1],A[2],A[3],A[4]),
                (x0,x1,A[2],y0,A[4]),(x0,x1,y1,A[3],A[4])]
    return [r for r in candidates if r[0]<r[1] and r[2]<r[3]]


def interval(lo,hi,D2,L):
    assert 0<=lo<=hi
    return {'lower':str(F(lo,D2)),'upper':str(F(hi,D2)),
            'normalized_lower':str(F(lo,D2)/L),'normalized_upper':str(F(hi,D2)/L)}


def evaluate(mod,case,block):
    p=[tuple(map(F,q)) for q in case['points']];m=list(map(F,case['masses']))
    S=block['S_mask'];L=F(block['L']);theta=F(block['theta']);bits=block['bits'];D=1<<bits;D2=D*D
    ids=[i for i in range(len(p)) if S>>i&1]
    a=min(p[i][0] for i in ids);b=max(p[i][0] for i in ids);c=min(p[i][1] for i in ids);d=max(p[i][1] for i in ids)
    coords=[q*D for q in (a,b,c,d)];assert all(q.denominator==1 for q in coords)
    aD,bD,cD,dD=map(int,coords);ri,ro=radii(theta*L,bits)
    ci=(aD-ri,bD+ri,cD-ri,dD+ri,S);co=(aD-ro,bD+ro,cD-ro,dD+ro,S)
    rings_i=[];rings_o=[];rows=[]
    for B in mod.closures(p):
        mask=B['mask'];mass=sum(m[i] for i in B['ids'])
        if not mask&S or mass<=theta*L or mass<=B['span2']:continue
        aa=sum(m[i] for i in B['ids'] if S>>i&1);bb=mass-aa
        ei,eo=receiver(B['bbox'],mass,bits,mask);bi,bo=receiver(B['bbox'],bb,bits,mask)
        ring_i=[q for r in ei for q in subtract_rect(r,bo)]
        ring_o=[q for r in eo for q in subtract_rect(r,bi)]
        rings_i+=ring_i;rings_o+=ring_o
        ar_i=mod.difference_area(ring_i,[]);ar_o=mod.difference_area(ring_o,[])
        outside_i=mod.difference_area(ring_i,[co]);outside_o=mod.difference_area(ring_o,[ci])
        assert F(ar_i,D2)<=aa # analytic ring cap; outward upper may exceed a by rounding.
        rows.append({'mask':mask,'rank':mask.bit_count(),'is_S':mask==S,'a_intersection_mass':str(aa),
                     'b_remaining_mass':str(bb),'total_mass':str(mass),
                     'erosion_ring_area':interval(ar_i,ar_o,D2,L),
                     'erosion_ring_area_analytic_upper_a':str(aa),
                     'erosion_ring_outside_container':interval(outside_i,outside_o,D2,L),
                     'geometric_ring_upper_at_most_a':F(ar_o,D2)<=aa,
                     'scope':'Individual bound only; mass intersections overlap between competing T.'})
    lo=mod.difference_area(rings_i,[co]);hi=mod.difference_area(rings_o,[ci])
    assert len(rows)==block['high_competitor_count_including_S']
    assert F(block['high_unpaid_union']['lower'])<=F(hi,D2)
    return {'S_mask':S,'L':str(L),'ring_union_outside_container':interval(lo,hi,D2,L),
            'actual_high_unpaid_union_from_budget':block['high_unpaid_union'],
            'competitors':rows,'ring_rectangle_pieces_inner':len(rings_i),'ring_rectangle_pieces_outer':len(rings_o),
            'scope':'Union of all exact-enclosed rings; upper bounds true high unpaid outside container; no sum(a) budget assertion.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',type=Path,required=True)
    ap.add_argument('--budget-result',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--case-id');args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert sha(args.module)==GEOMETRY_SHA
    spec=importlib.util.spec_from_file_location('frozen_r064_geometry',args.module)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    data=json.loads(args.budget_result.read_text());assert data['status']=='PASS' and data['source_sha256']==BUDGET_SHA
    started=time.perf_counter();result={'schema':'r066-erosion-ring-diagnostic-v1','status':'RUNNING',
                                       'source_sha256':sha(Path(__file__)),'geometry_sha256':sha(args.module),
                                       'budget_result_sha256':sha(args.budget_result),'python':sys.version,'cases':[]}
    for case in data['cases']:
        if args.case_id and case['id']!=args.case_id:continue
        result['cases'].append({'id':case['id'],'blocks':[evaluate(mod,case,b) for b in case['blocks']]})
    result['status']='PASS';result['wall_seconds']=time.perf_counter()-started
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cases':len(result['cases']),'seconds':result['wall_seconds']}))


if __name__=='__main__':main()
