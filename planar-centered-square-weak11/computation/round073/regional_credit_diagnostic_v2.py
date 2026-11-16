#!/usr/bin/env python3
"""R73 exact clipped-region bookkeeping from frozen R72 halo enclosures.

No source closure, square root or old complete halo is recomputed. New
vertical x-event interval-merging sweeps clip the certified inner/outer rectangles
to a disjoint 3 by 3 C/J partition. Unit-window masses are enumerated exactly
over all finite source-coordinate origin strata, including boundary points.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,sys,time
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

CASES=[('s25-only','s25-only.json'),('r071-positive33','r071-positive33.json'),('aligned225','aligned225-actual-input.json'),('contracted-unique','contracted-unique-actual-input.json')]
CL=F(-1,2);CH=F(73,50);JL=F(-1,25);JH=F(1)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def merge(iv):
    out=[]
    for a,b in sorted(iv):
      if not out or a>out[-1][1]:out.append([a,b])
      elif b>out[-1][1]:out[-1][1]=b
    return out
def intersection_length(a,b):
    i=j=0;out=0
    while i<len(a) and j<len(b):
      lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
      if lo<hi:out+=hi-lo
      if a[i][1]<=b[j][1]:i+=1
      else:j+=1
    return out
def sweep(left,right):
    """Integer area of union(left), union(right), and left-minus-right."""
    events=defaultdict(list)
    for channel,rr in enumerate((left,right)):
      for j,(a,b,c,d) in enumerate(rr):
        if a>=b or c>=d:continue
        events[a].append((channel,j,(c,d)));events[b].append((channel,j,None))
    if not events:return (0,0,0)
    active=[{},{}];area=[0,0,0];last=min(events)
    for x,ee in sorted(events.items()):
      if x>last:
        a=merge(active[0].values());b=merge(active[1].values());la=sum(d-c for c,d in a);lb=sum(d-c for c,d in b)
        area[0]+=(x-last)*la;area[1]+=(x-last)*lb;area[2]+=(x-last)*(la-intersection_length(a,b))
      for channel,j,item in ee:
        if item is None:del active[channel][j]
        else:active[channel][j]=item
      last=x
    assert not active[0] and not active[1]
    return tuple(area)
def clip(rr,box):
    u,v,s,t=box;out=[]
    for a,b,c,d in rr:
      z=(max(a,u),min(b,v),max(c,s),min(d,t))
      if z[0]<z[1] and z[2]<z[3]:out.append(z)
    return out
def enc(lo,hi):
    lo=F(lo);hi=F(hi);assert 0<=lo<=hi
    return {'lower':str(lo),'upper':str(hi),'lower_decimal':float(lo),'upper_decimal':float(hi),'width':str(hi-lo)}
def aggregate(records,ids):
    keys=['full_halo','residual_halo','deletion','missing_deletion','full_halo_vacancy'];out={'region_ids':ids,'area':str(sum(F(records[k]['area']) for k in ids))}
    for metric in keys:out[metric]=enc(sum(F(records[k][metric]['lower']) for k in ids),sum(F(records[k][metric]['upper']) for k in ids))
    return out

def window_diagnostic(data):
    points=[tuple(map(F,p)) for p in data['positions']];masses=list(map(F,data['masses']));S=set(data['S_indices']);sx=sorted(set(points[k][0] for k in S));sy=sorted(set(points[k][1] for k in S));assert len(sx)==len(sy)==5
    axes=[]
    for axis in (0,1):
      edges=sorted({p[axis] for p in points}|{p[axis]-1 for p in points});reps=[]
      for i,a in enumerate(edges):
        reps.append((a,'point'))
        if i+1<len(edges):reps.append(((a+edges[i+1])/2,'open_interval'))
      axes.append(reps)
    groups={};families={};masks={};sat={};strata=0;window_count=0
    def summary(container,key,row):
      if key not in container or F(row['slack'])<F(container[key]['slack']):container[key]=row
    def kind(captured,grid,low,high):
      if len(captured)==5:return 'all'
      if captured and captured[0]==grid[0]:return low
      if captured and captured[-1]==grid[-1]:return high
      return 'empty'
    for (x,xt) in axes[0]:
      xc=[a for a in sx if x<=a<=x+1];xk=kind(xc,sx,'left','right')
      for y,yt in axes[1]:
        strata+=1;yc=[a for a in sy if y<=a<=y+1];yk=kind(yc,sy,'bottom','top');ids=[i for i,p in enumerate(points) if x<=p[0]<=x+1 and y<=p[1]<=y+1]
        if not ids:continue
        window_count+=1;mask=sum(1<<i for i in ids);mass=sum(masses[i] for i in ids);inside=sum(masses[i] for i in ids if i in S);outside=mass-inside;k=len(xc);l=len(yc);assert inside==F(k*l,25) and mass<=1
        row={'k':k,'l':l,'internal_mass':str(inside),'outside_mass':str(outside),'outside_capacity':str(1-inside),'total_mass':str(mass),'slack':str(1-mass),'capture_mask':str(mask),'captured_indices':ids,'representative_closed_unit_square':[str(x),str(x+1),str(y),str(y+1)],'origin_stratum_types':[xt,yt]}
        masks.setdefault(mask,row)
        summary(groups,f'{k}x{l}',row)
        if xk=='all' and yk=='all':family='all_S'
        elif xk in ('left','right') and yk=='all':family='edge_'+xk
        elif yk in ('bottom','top') and xk=='all':family='edge_'+yk
        elif xk in ('left','right') and yk in ('bottom','top'):family='corner_'+xk+'_'+yk
        else:family='no_internal_'+xk+'_'+yk
        summary(families,family,row)
        if mass==1:sat.setdefault(mask,row)
    return {'method':'All point/open-interval strata of unit-window origins at source_coordinate and source_coordinate-1; every nonempty finite capture pattern represented, including closed boundary coincidences.','axis_origin_strata_counts':list(map(len,axes)),'cartesian_origin_strata':strata,'nonempty_origin_strata':window_count,'distinct_nonempty_capture_masks':len(masks),'saturating_unit_window_capture_mask_count':len(sat),'saturating_unit_window_masks':list(sat.values()),'minimum_slack_by_internal_kxl':groups,'minimum_slack_by_side_corner_family':families,'global_minimum_unit_window_slack':str(min(F(r['slack']) for r in masks.values()))}

def run_case(root,name,input_name):
    start=time.monotonic();gp=root/'computation/results/round072-actual'/f'{name}.json';ip=root/'computation/inputs/r072'/input_name;d=json.loads(gp.read_text());data=json.loads(ip.read_text());assert sha(ip)==d['receipt']['input_sha256'];job=next(a for a in d['all_maximum_block_evaluations'] if a['chosen_S']);D=d['coordinate_scale'];den=D*D
    fi=[tuple(r['inner']) for r in d['full_directed_rectangles']];fo=[tuple(r['outer']) for r in d['full_directed_rectangles']];ri=[tuple(r['inner']) for r in job['residual_directed_rectangles']];ro=[tuple(r['outer']) for r in job['residual_directed_rectangles']]
    intervals=[(CL,JL),(JL,JH),(JH,CH)];records={};corners=[];sides=[]
    for i in range(3):
      for j in range(3):
        tag=f'x{i}y{j}';a,b=intervals[i];c,e=intervals[j];bb=[a,b,c,e];box=[]
        for v in bb:assert (v*D).denominator==1;box.append(int(v*D))
        area=(b-a)*(e-c);lf,hr,ld=sweep(clip(fi,box),clip(ro,box));hf,lr,hd=sweep(clip(fo,box),clip(ri,box));assert 0<=ld<=hd<=area*den
        record={'bounds':[str(v) for v in bb],'axis_interval_types':[('left_closed_right_open','closed','left_open_right_closed')[i],('left_closed_right_open','closed','left_open_right_closed')[j]],'area':str(area),'kind':'J' if i==j==1 else 'corner' if i!=1 and j!=1 else 'side','full_halo':enc(F(lf,den),F(hf,den)),'residual_halo':enc(F(lr,den),F(hr,den)),'deletion':enc(F(ld,den),F(hd,den)),'missing_deletion':enc(area-F(hd,den),area-F(ld,den)),'full_halo_vacancy':enc(area-F(hf,den),area-F(lf,den)),'clipped_rectangle_counts':{'full_inner':len(clip(fi,box)),'full_outer':len(clip(fo,box)),'residual_inner':len(clip(ri,box)),'residual_outer':len(clip(ro,box))}}
        # True missing = vacancy plus residual, since H(nu) subset H(mu).
        lhs=(F(record['missing_deletion']['lower']),F(record['missing_deletion']['upper']));rhs=(F(record['full_halo_vacancy']['lower'])+F(record['residual_halo']['lower']),F(record['full_halo_vacancy']['upper'])+F(record['residual_halo']['upper']));assert max(lhs[0],rhs[0])<=min(lhs[1],rhs[1])
        records[tag]=record
        if record['kind']=='corner':corners.append(tag)
        elif record['kind']=='side':sides.append(tag)
    aggregates={'C':aggregate(records,list(records)),'J':aggregate(records,['x1y1']),'C_minus_J':aggregate(records,[k for k in records if k!='x1y1']),'corners':aggregate(records,corners),'sides':aggregate(records,sides)}
    old=(F(job['actual_Delta']['lower']),F(job['actual_Delta']['upper']));new=(F(aggregates['C']['deletion']['lower']),F(aggregates['C']['deletion']['upper']));assert new[0]>=old[0] and new[1]<=old[1] # Clipping may tighten only roundoff bands outside exact container C.
    target=(CH-CL)**2-F(7,2);cornerthreshold=F(aggregates['corners']['area'])-target
    diagnostics={'C_area':str((CH-CL)**2),'J_area':str((JH-JL)**2),'target_missing_credit_Carea_minus_7over2':str(target),'outer_ring_missing_lower_minus_target':str(F(aggregates['C_minus_J']['missing_deletion']['lower'])-target),'outer_ring_candidate_proved_for_this_control':F(aggregates['C_minus_J']['missing_deletion']['lower'])>=target,'corner_deletion_threshold_area_minus_target':str(cornerthreshold),'corner_missing_lower_minus_target':str(F(aggregates['corners']['missing_deletion']['lower'])-target),'corner_only_candidate_proved_for_this_control':F(aggregates['corners']['missing_deletion']['lower'])>=target,'missing_credit_decomposition':'For each disjoint region, missing deletion = full-halo vacancy + residual halo. Every quantity is a directed rational enclosure, not an exact radical evaluation.'}
    windows=window_diagnostic(data)
    return {'name':name,'input_path':str(ip),'input_sha256':sha(ip),'certified_geometry_path':str(gp),'certified_geometry_sha256':sha(gp),'chosen_S_mask':job['maximum_mask'],'source_geometry_bits':d['bits'],'coordinate_scale':D,'positive_source_count':d['positive_source_count'],'maximum_count':d['maximum_count'],'nonterminal':d['nonterminal'],'regions':records,'aggregates':aggregates,'diagnostics':diagnostics,'unit_window_caps':windows,'seconds':time.monotonic()-start}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();t=time.monotonic();root=Path(a.archive);cases=[run_case(root,*p) for p in CASES]
    out={'schema':'r073-clipped-region-and-unit-cap-diagnostic-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'source_archive':str(root),'parameters':{'C':['-1/2','73/50','-1/2','73/50'],'J':['-1/25','1','-1/25','1'],'partition':'three disjoint axis intervals [Clo,Jlo),[Jlo,Jhi],(Jhi,Chi]; their9 Cartesian products partition closedC exactly. Area sweeps may use closed bounds because boundaries have measure0.','random_seed':None,'new_root_computations':0,'new_closure_enumerations':0},'cases':cases,'seconds':time.monotonic()-t,'scope':'Four fixed actual atomic controls only. New exact clipped-region sweeps reuse certified72bit full/residual rectangles. Unit-window source-coordinate strata are complete finite exact mass diagnostics. No general corner inequality, new actual input, MIP or general C2 bound is inferred.'}
    Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'output':a.out,'pid':os.getpid(),'seconds':out['seconds'],'cases':[{'name':c['name'],'Delta_C':c['aggregates']['C']['deletion'],'corner_Delta':c['aggregates']['corners']['deletion'],'corner_credit_margin':c['diagnostics']['corner_missing_lower_minus_target'],'unit_saturating_masks':c['unit_window_caps']['saturating_unit_window_capture_mask_count']} for c in cases]},indent=2))

if __name__=='__main__':main()
