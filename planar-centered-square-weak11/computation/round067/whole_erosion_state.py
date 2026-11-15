#!/usr/bin/env python3
r"""Mass-sensitive full erosion-ring unions W and Wprime, with no theta container.

All strictly eligible full-source bbox closures intersecting S are retained.
W uses the full-T bbox at outside mass b; Wprime uses the tight bbox of T\S.
Neither state's area is discounted by the full remaining halo. Hremaining is
computed separately only to diagnose true deletion and overpayment.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def bounds(lo,hi,D2,L):
    assert 0<=lo<=hi
    return {'lower':str(F(lo,D2)),'upper':str(F(hi,D2)),
            'normalized_lower':str(F(lo,D2)/L),'normalized_upper':str(F(hi,D2)/L)}


def intersection_area(mod,left,right):
    return mod.difference_area(left,[])-mod.difference_area(left,right)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--module',type=Path,required=True);ap.add_argument('--erosion-source',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--bits',type=int,default=72)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data') and args.bits>=16
    assert sha(args.module)=='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
    assert sha(args.erosion_source)=='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
    started=time.perf_counter();mod=load('frozen_r064_geometry',args.module);er=load('frozen_r066_rectclip',args.erosion_source)
    data=json.loads(args.input.read_text());p=[tuple(map(F,q)) for q in data['points']]
    mass=list(map(F,data['initial_masses']));n=len(p);S=data['S_mask'];assert 0<S<(1<<n)-1
    assert len(mass)==n and min(mass)>0 and len(set(p))==n
    L=sum(mass[i] for i in range(n) if S>>i&1);D2=1<<(2*args.bits)
    boxes=mod.closures(p);eligible=[b for b in boxes if sum(mass[i] for i in b['ids'])>b['span2']]
    maximum=max(sum(mass[i] for i in b['ids']) for b in eligible)
    maximum_masks=[b['mask'] for b in eligible if sum(mass[i] for i in b['ids'])==maximum]
    assert L==maximum and S in maximum_masks
    full_i,full_o=mod.exact_rectangles(p,mass,boxes,args.bits)
    fi={r[4]:[r] for r in full_i};fo={r[4]:[r] for r in full_o}
    wi=[];wo=[];pi=[];po=[];rows=[]
    for box in eligible:
        mask=box['mask']
        if not mask&S:continue
        total=sum(mass[i] for i in box['ids']);outside=[i for i in box['ids'] if not S>>i&1]
        b=sum((mass[i] for i in outside),F(0));a=total-b
        ki,ko=er.receiver(box['bbox'],b,args.bits,mask)
        if outside:
            tight=(min(p[i][0] for i in outside),max(p[i][0] for i in outside),
                   min(p[i][1] for i in outside),max(p[i][1] for i in outside))
            ti,to=er.receiver(tight,b,args.bits,mask)
            tight_span2=max(tight[1]-tight[0],tight[3]-tight[2])**2
        else:ti,to=[],[];tight_span2=None
        ii=fi.get(mask,[]);oo=fo.get(mask,[])
        rwi=[q for rect in ii for q in er.subtract_rect(rect,ko)]
        rwo=[q for rect in oo for q in er.subtract_rect(rect,ki)]
        rpi=[q for rect in ii for q in er.subtract_rect(rect,to)]
        rpo=[q for rect in oo for q in er.subtract_rect(rect,ti)]
        assert mod.difference_area(rpi,rwi)==0 and mod.difference_area(rpo,rwo)==0
        wi+=rwi;wo+=rwo;pi+=rpi;po+=rpo
        wlo=mod.difference_area(rwi,[]);whi=mod.difference_area(rwo,[])
        plo=mod.difference_area(rpi,[]);phi=mod.difference_area(rpo,[])
        assert F(wlo,D2)<=a and F(plo,D2)<=a
        rows.append({'mask':mask,'rank':mask.bit_count(),'intersect_S_rank':(mask&S).bit_count(),
                     'mass':str(total),'a_intersection_mass':str(a),'b_outside_mass':str(b),
                     'full_bbox_span_squared':str(box['span2']),
                     'outside_tight_bbox_span_squared':str(tight_span2) if tight_span2 is not None else None,
                     'full_bbox_kernel_nonempty':box['span2']<b,
                     'tight_outside_kernel_nonempty':bool(outside) and tight_span2<b,
                     'W_label_ring':bounds(wlo,whi,D2,L),'Wprime_label_ring':bounds(plo,phi,D2,L),
                     'analytic_single_label_upper_a':str(a)})
    keep=[i for i in range(n) if not S>>i&1];rp=[p[i] for i in keep];rm=[mass[i] for i in keep]
    rem_i,rem_o=mod.exact_rectangles(rp,rm,mod.closures(rp),args.bits)
    wlo=mod.difference_area(wi,[]);whi=mod.difference_area(wo,[])
    plo=mod.difference_area(pi,[]);phi=mod.difference_area(po,[])
    delta_lo=mod.difference_area(full_i,rem_o);delta_hi=mod.difference_area(full_o,rem_i)
    # Finite integer inclusion checks, stronger than interval overlap alone.
    checks={'Wprime_inner_minus_W_inner':mod.difference_area(pi,wi),
            'Wprime_outer_minus_W_outer':mod.difference_area(po,wo),
            'W_outer_minus_Hfull_outer':mod.difference_area(wo,full_o),
            'Hfull_inner_minus_union_Hrem_outer_Wprime_inner':mod.difference_area(full_i,rem_o+pi)}
    assert all(v==0 for v in checks.values()),checks
    # Since Delta subset Wprime subset W subset Hfull, W-Delta=W intersect
    # Hremaining. These intersections diagnose overpayment; they DO NOT alter
    # the W/Wprime state costs above.
    owlo=intersection_area(mod,wi,rem_i);owhi=intersection_area(mod,wo,rem_o)
    oplo=intersection_area(mod,pi,rem_i);ophi=intersection_area(mod,po,rem_o)
    residual_w=(wlo-delta_hi-owhi,whi-delta_lo-owlo)
    residual_p=(plo-delta_hi-ophi,phi-delta_lo-oplo)
    assert residual_w[0]<=0<=residual_w[1] and residual_p[0]<=0<=residual_p[1]
    refinement_lo=mod.difference_area(wi,po);refinement_hi=mod.difference_area(wo,pi)
    payload={'schema':'r067-whole-erosion-state-v1','status':'PASS','id':data['id'],
             'python':sys.version,'source_sha256':sha(Path(__file__)),'input_sha256':sha(args.input),
             'geometry_source_sha256':sha(args.module),'rectclip_source_sha256':sha(args.erosion_source),
             'bits':args.bits,'source_count':n,'S_mask':S,'S_rank':S.bit_count(),'L':str(L),
             'closed_bbox_count':len(boxes),'strict_eligible_count':len(eligible),'intersect_S_eligible_count':len(rows),
             'maximum_eligible_mass':str(maximum),'all_maximum_masks':maximum_masks,
             'W':bounds(wlo,whi,D2,L),'Wprime':bounds(plo,phi,D2,L),
             'true_deletion':bounds(delta_lo,delta_hi,D2,L),
             'full_halo':bounds(mod.difference_area(full_i,[]),mod.difference_area(full_o,[]),D2,L),
             'remaining_halo':bounds(mod.difference_area(rem_i,[]),mod.difference_area(rem_o,[]),D2,L),
             'W_minus_Delta_overpayment':bounds(owlo,owhi,D2,L),
             'Wprime_minus_Delta_overpayment':bounds(oplo,ophi,D2,L),
             'W_minus_Wprime_refinement':bounds(refinement_lo,refinement_hi,D2,L),
             'overpayment_identity_W_residual':list(map(lambda v:str(F(v,D2)),residual_w)),
             'overpayment_identity_Wprime_residual':list(map(lambda v:str(F(v,D2)),residual_p)),
             'exact_zero_area_inclusion_checks':{k:str(v) for k,v in checks.items()},
             'certifies_W_le_7over2_L':F(whi,D2)<=F(7,2)*L,'certifies_Wprime_le_7over2_L':F(phi,D2)<=F(7,2)*L,
             'all_label_rings':rows,'ring_piece_counts':{'W_inner':len(wi),'W_outer':len(wo),'Wprime_inner':len(pi),'Wprime_outer':len(po)},
             'wall_seconds':time.perf_counter()-started,
             'scope':'Two mass-sensitive state costs use all eligible closures intersecting S; no theta, no complete-Hremaining discount in their construction. Remaining halo is used only for separate deletion/overpayment diagnostics. Input-specific certificates only.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'status':'PASS','id':data['id'],'labels':len(rows),'W_upper_display':float(F(whi,D2)/L),
                      'Wprime_upper_display':float(F(phi,D2)/L),'Delta_upper_display':float(F(delta_hi,D2)/L)}))


if __name__=='__main__':main()
