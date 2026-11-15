#!/usr/bin/env python3
"""Exact low-container/high-competitor union budget on prescribed finite inputs.

No LP and no random samples. Root owns execution. High competitors include S:
E_S need not lie in C_theta for arbitrary source geometry.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import sys
import time


MODULE_SHA='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def root_radius(mass,bits):
    """Strict dyadic radii below/above sqrt(mass)/2, exact squared checks."""
    D=1<<bits;value=mass.numerator*D*D;den=4*mass.denominator
    k=math.isqrt(value//den)
    assert den*k*k<=value<den*(k+1)*(k+1)
    return max(0,k-1),k+1


def bbox(points,ids):
    return (min(points[i][0] for i in ids),max(points[i][0] for i in ids),
            min(points[i][1] for i in ids),max(points[i][1] for i in ids))


def area(rect):return max(0,rect[1]-rect[0])*max(0,rect[3]-rect[2])


def interval(lo,hi,D2,L=F(1)):
    assert 0<=lo<=hi,(lo,hi)
    return {'lower':str(F(lo,D2)),'upper':str(F(hi,D2)),
            'normalized_lower':str(F(lo,D2)/L),'normalized_upper':str(F(hi,D2)/L)}


def components(points,masses,eligible):
    cs=[1<<i for i in range(len(points))]
    for mask in eligible:
        joined=mask
        for c in cs:
            if c&mask:joined |= c
        cs=[c for c in cs if not c&mask]+[joined]
    result=[]
    for mask in sorted(cs):
        ids=[i for i in range(len(points)) if mask>>i&1];a,b,c,d=bbox(points,ids)
        mass=sum(masses[i] for i in ids);span2=max(b-a,d-c)**2
        result.append({'mask':mask,'rank':mask.bit_count(),'mass':str(mass),
                       'span_squared_over_mass':str(span2/mass),'terminal_35_33':33*span2<=35*mass})
    return result


def evaluate(mod,points,masses,S,theta,bits,with_individuals=True):
    started=time.perf_counter();n=len(points);assert 0<S<(1<<n)-1
    ids=[i for i in range(n) if S>>i&1];L=sum(masses[i] for i in ids)
    boxes=mod.closures(points);massmap={b['mask']:sum(masses[i] for i in b['ids']) for b in boxes}
    eligible=[b for b in boxes if massmap[b['mask']]>b['span2']]
    maximum=max(massmap[b['mask']] for b in eligible)
    assert L==maximum and S in {b['mask'] for b in eligible}
    high=[b for b in eligible if b['mask']&S and massmap[b['mask']]>theta*L]
    high_masks={b['mask'] for b in high};assert S in high_masks
    full_i,full_o=mod.exact_rectangles(points,masses,boxes,bits)
    hi_i=[r for r in full_i if r[4] in high_masks];hi_o=[r for r in full_o if r[4] in high_masks]
    keep=[i for i in range(n) if not S>>i&1]
    rp=[points[i] for i in keep];rm=[masses[i] for i in keep]
    rem_i,rem_o=mod.exact_rectangles(rp,rm,mod.closures(rp),bits)
    D=1<<bits;D2=D*D;a,b,c,d=bbox(points,ids)
    coords=[q*D for q in (a,b,c,d)];assert all(q.denominator==1 for q in coords)
    aD,bD,cD,dD=map(int,coords);ri,ro=root_radius(theta*L,bits)
    ci=(aD-ri,bD+ri,cD-ri,dD+ri,S);co=(aD-ro,bD+ro,cD-ro,dD+ro,S)
    ca_i,ca_o=area(ci),area(co)
    old_i=ca_i-mod.difference_area([ci],rem_i)
    old_o=ca_o-mod.difference_area([co],rem_o)
    high_i=mod.difference_area(hi_i,[co]+rem_o)
    high_o=mod.difference_area(hi_o,[ci]+rem_i)
    raw_hi_i=mod.difference_area(hi_i,[co]);raw_hi_o=mod.difference_area(hi_o,[ci])
    hole_i=mod.difference_area([ci],full_o);hole_o=mod.difference_area([co],full_i)
    delta_i=mod.difference_area(full_i,rem_o);delta_o=mod.difference_area(full_o,rem_i)
    budget_i=max(0,ca_i-old_o+high_i);budget_o=ca_o-old_i+high_o
    # The analytic identity is core's theorem. These are necessary exact
    # consistency checks between independent monotone enclosure expressions.
    residual_i=budget_i-delta_o-hole_o;residual_o=budget_o-delta_i-hole_i
    assert residual_i<=0<=residual_o,(residual_i,residual_o)
    assert budget_o>=delta_i
    output={'S_mask':S,'S_rank':len(ids),'L':str(L),'theta':str(theta),'bits':bits,
            'closed_bbox_count':len(boxes),'strict_eligible_count':len(eligible),
            'maximum_eligible_mass':str(maximum),'high_competitor_count_including_S':len(high),
            'S_bbox':list(map(str,(a,b,c,d))),'S_width_squared_over_L':[str((b-a)**2/L),str((d-c)**2/L)],
            'container_area':interval(ca_i,ca_o,D2,L),'low_old_credit':interval(old_i,old_o,D2,L),
            'high_unpaid_union':interval(high_i,high_o,D2,L),
            'high_union_outside_container_before_old_credit':interval(raw_hi_i,raw_hi_o,D2,L),
            'low_hole':interval(hole_i,hole_o,D2,L),'true_deletion':interval(delta_i,delta_o,D2,L),
            'analytic_envelope_budget':interval(budget_i,budget_o,D2,L),
            'identity_budget_minus_delta_minus_hole':{'lower':str(F(residual_i,D2)),'upper':str(F(residual_o,D2)),'contains_zero':True},
            'high_unpaid_minus_low_old_credit_normalized':{'lower':str(F(high_i-old_o,D2)/L),'upper':str(F(high_o-old_i,D2)/L)},
            'budget_certifies_7over2':F(budget_o,D2)<=F(7,2)*L,
            'components':components(points,masses,[b['mask'] for b in eligible]),
            'scope':'Exact input-specific intervals. High geometry aggregated as a union, never as a sum of individual areas.'}
    if theta==F(3,4):
        # (7/4 - sqrt3) is positive; rational outward enclosure from integer
        # floor sqrt(3)*2^bits. Compare with its LOWER endpoint for sufficiency.
        k=math.isqrt(3*D2);assert k*k<=3*D2<(k+1)*(k+1)
        threshold_lo=F(7,4)-F(k+1,D);threshold_hi=F(7,4)-F(k,D)
        output['universal_compensation_threshold']={'lower':str(threshold_lo),'upper':str(threshold_hi)}
        output['certifies_universal_compensation_condition']=F(high_o-old_i,D2)/L<=threshold_lo
    if with_individuals:
        imap={r[4]:r for r in hi_i};omap={r[4]:r for r in hi_o};rows=[]
        for binfo in high:
            mask=binfo['mask'];ir=[imap[mask]] if mask in imap else [];orr=[omap[mask]] if mask in omap else []
            ar_i=sum(area(r) for r in ir);ar_o=sum(area(r) for r in orr)
            u_i=mod.difference_area(ir,[co]+rem_o);u_o=mod.difference_area(orr,[ci]+rem_i)
            outside_i=mod.difference_area(ir,[co]);outside_o=mod.difference_area(orr,[ci])
            aa,bb,cc,dd=binfo['bbox'];m=massmap[mask]
            rows.append({'mask':mask,'rank':mask.bit_count(),'is_S':mask==S,'mass_over_L':str(m/L),
                         'bbox_span_squared_over_mass':str(binfo['span2']/m),
                         'bbox_width_squared_over_mass':[str((bb-aa)**2/m),str((dd-cc)**2/m)],
                         'receiver_area':interval(ar_i,ar_o,D2,L),
                         'outside_container_area':interval(outside_i,outside_o,D2,L),
                         'individually_unpaid_after_remaining_halo':interval(u_i,u_o,D2,L),
                         'certified_contained_in_container':outside_o==0,
                         'certified_individually_fully_paid':u_o==0,
                         'certified_receiver_area_at_most_L_over_10000':F(ar_o,D2)<=L/10000})
        output['individual_competitors_DIAGNOSTIC_ONLY_NOT_UNION_SUM']=rows
    output['wall_seconds']=time.perf_counter()-started
    return output


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',type=Path,required=True)
    inp=ap.add_mutually_exclusive_group(required=True);inp.add_argument('--r64-input',type=Path);inp.add_argument('--r63-input',type=Path)
    ap.add_argument('--stratum',type=int,default=0);ap.add_argument('--case',type=int)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--bits',type=int,default=72)
    ap.add_argument('--theta',default='3/4');ap.add_argument('--max-blocks',type=int)
    ap.add_argument('--seconds-per-block',type=int,default=120);ap.add_argument('--total-seconds',type=int,default=900)
    ap.add_argument('--no-individuals',action='store_true');args=ap.parse_args()
    assert os.path.ismount('/Volumes/data') and sha(args.module)==MODULE_SHA
    assert args.bits>=16 and 0<F(args.theta)<1 and args.seconds_per_block>0 and args.total_seconds>0
    spec=importlib.util.spec_from_file_location('frozen_r064_geometry',args.module)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    source=args.r64_input or args.r63_input;data=json.loads(source.read_text());cases=[]
    if args.r64_input:
        cases=[{'id':data['id'],'points':[tuple(map(F,p)) for p in data['points']],
                'masses':list(map(F,data['initial_masses'])),'selected_S':[data['S_mask']]}]
    else:
        for original in data['strata'][args.stratum]['inputs']:
            p=[(F(x,4),F(y,4)) for x,y in original['points_quarters']];m=[F(v,16) for v in original['mass_sixteenths']]
            bs=mod.closures(p);eligible=[b for b in bs if sum(m[i] for i in b['ids'])>b['span2']]
            largest=max(sum(m[i] for i in b['ids']) for b in eligible)
            masks=[b['mask'] for b in eligible if sum(m[i] for i in b['ids'])==largest]
            cases.append({'id':original['id'],'points':p,'masses':m,'selected_S':masks})
    if args.case is not None:cases=[cases[args.case]]
    result={'schema':'r066-theta-competitor-budget-v1','status':'RUNNING','source_sha256':sha(Path(__file__)),
            'geometry_module_sha256':sha(args.module),'input_sha256':sha(source),'input_path':str(source),
            'python':sys.version,'theta':str(F(args.theta)),'cases':[],
            'budget':{'seconds_per_block':args.seconds_per_block,'total_seconds':args.total_seconds},
            'scope':'Frozen deterministic inputs and all proper maximum blocks. Exact integer/Fraction enclosures; no new universal theorem or C2 bound.'}
    def alarm(signum,frame):raise TimeoutError('block wall budget exceeded')
    signal.signal(signal.SIGALRM,alarm);started=time.monotonic()
    for case in cases:
        r={'id':case['id'],'points':[[str(x),str(y)] for x,y in case['points']],
           'masses':list(map(str,case['masses'])),'all_maximum_masks':case['selected_S'],'blocks':[]}
        masks=[s for s in case['selected_S'] if s!=(1<<len(case['points']))-1]
        r['improper_full_maximum_masks_skipped']=[s for s in case['selected_S'] if s not in masks]
        if args.max_blocks is not None:masks=masks[:args.max_blocks]
        r['all_proper_maximum_masks_evaluated']=len(masks)==len(case['selected_S'])-len(r['improper_full_maximum_masks_skipped'])
        result['cases'].append(r)
        try:
            for S in masks:
                remaining=args.total_seconds-(time.monotonic()-started)
                if remaining<=0:raise TimeoutError('total wall budget exceeded')
                signal.alarm(max(1,min(args.seconds_per_block,int(remaining))))
                r['blocks'].append(evaluate(mod,case['points'],case['masses'],S,F(args.theta),args.bits,not args.no_individuals))
                signal.alarm(0);write(args.out,result)
            write(args.out,result)
        except Exception as exc:
            signal.alarm(0);result['status']='FAILED';result['error']=repr(exc);write(args.out,result);raise
    result['status']='PASS';result['wall_seconds']=time.monotonic()-started
    result['evaluated_proper_blocks']=sum(len(r['blocks']) for r in result['cases']);write(args.out,result)
    print(json.dumps({'status':'PASS','cases':len(result['cases']),'blocks':result['evaluated_proper_blocks'],'seconds':result['wall_seconds']}))


if __name__=='__main__':main()
