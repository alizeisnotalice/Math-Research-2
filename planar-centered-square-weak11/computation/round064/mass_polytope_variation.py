#!/usr/bin/env python3
"""One fixed-position LP proposal, followed by exact full-halo certification.

Finite differences and SciPy LP are numerical proposals only, including at the
symmetric seed's walls. Every accepted mass vector and every claimed positive
improvement is certified with Fraction/integers on complete bbox geometry.
No R63 code is imported or modified.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def closures(points):
    """Complete tight bbox source closures by coordinate-interval masks."""
    axes=[]
    for a in (0,1):
        groups=defaultdict(int)
        for i,p in enumerate(points):groups[p[a]] |= 1<<i
        xs=sorted(groups);masks=[]
        for j in range(len(xs)):
            mask=0
            for k in range(j,len(xs)):mask |= groups[xs[k]];masks.append(mask)
        axes.append(masks)
    result=[]
    for mask in sorted({a&b for a in axes[0] for b in axes[1]}-{0}):
        ids=[i for i in range(len(points)) if mask>>i&1]
        bounds=[min(points[i][a] for i in ids) for a in (0,1)]+[max(points[i][a] for i in ids) for a in (0,1)]
        a,c,b,d=bounds;span=max(b-a,d-c)
        result.append({'mask':mask,'ids':ids,'bbox':(a,b,c,d),'span2':span*span})
    return result


def constraints(spec,boxes):
    n=len(spec['points']);S=spec['S_mask'];eta=F(spec['parameters']['unique_margin_eta'])
    C=spec['mass_constraints'];rows=[]
    for b in boxes:
        cap=max(F(1),b['span2'])
        rows.append((b['mask'],cap,'full-bbox-cap'))
        if b['mask']!=S and b['span2']<1:
            rows.append((b['mask'],1-eta,'unique-max-gap'))
    rows.append(((1<<n)-1,F(C['nonterminal_total_mass_cap']),'nonterminal-total-cap'))
    for edge in C['prescribed_connection_lower_constraints']:
        rows.append((edge['mask'],-F(edge['minimum_mass']),'negative-connection-lower'))
    return rows


def exact_check(masses,spec,rows,boxes):
    C=spec['mass_constraints'];floor=F(C['positive_floor_for_initial_chamber'])
    errors=[]
    if sum(masses[i] for i in spec['S_indices'])!=F(C['S_sum']):errors.append('S equality')
    if any(m<floor for m in masses):errors.append('positive floor')
    for mask,rhs,kind in rows:
        value=sum(masses[i] for i in range(len(masses)) if mask>>i&1)
        if kind=='negative-connection-lower':value=-value
        if value>rhs:errors.append(f'{kind}:{mask}')
    eligible=[b for b in boxes if sum(masses[i] for i in b['ids'])>b['span2']]
    top=max(sum(masses[i] for i in b['ids']) for b in eligible)
    winners=[b['mask'] for b in eligible if sum(masses[i] for i in b['ids'])==top]
    if top!=1 or winners!=[spec['S_mask']]:errors.append('maximum eligible mass/mask')
    return {'feasible':not errors,'errors':errors,'strict_eligible_closures':len(eligible),
            'maximum_eligible_mass':str(top),'maximum_eligible_masks':winners}


def merged(intervals):
    answer=[]
    for a,b in sorted(intervals):
        if answer and a<=answer[-1][1]:answer[-1]=(answer[-1][0],max(answer[-1][1],b))
        else:answer.append((a,b))
    return answer


def difference_area(left,right):
    ys=sorted({r[j] for r in left+right for j in (2,3)})
    total=0
    for a,b in zip(ys,ys[1:]):
        l=merged([(r[0],r[1]) for r in left if r[2]<=a and b<=r[3]])
        r=merged([(r[0],r[1]) for r in right if r[2]<=a and b<=r[3]])
        length=sum(d-c for c,d in l);i=j=0
        while i<len(l) and j<len(r):
            length-=max(0,min(l[i][1],r[j][1])-max(l[i][0],r[j][0]))
            if l[i][1]<r[j][1]:i+=1
            else:j+=1
        total+=(b-a)*length
    return total


def exact_rectangles(points,masses,boxes,bits):
    D=1<<bits;inner=[];outer=[]
    for b in boxes:
        mass=sum(masses[i] for i in b['ids'])
        if mass<=b['span2']:continue
        a,c,e,f=b['bbox'];coords=[v*D for v in (a,c,e,f)]
        assert all(v.denominator==1 for v in coords),'bits must cover dyadic source coordinates'
        a,c,e,f=map(int,coords)
        value=mass.numerator*D*D;den=4*mass.denominator
        k=math.isqrt(value//den)
        assert den*k*k<=value<den*(k+1)*(k+1)
        for rs,r in ((inner,max(0,k-1)),(outer,k+1)):
            rect=(c-r,a+r,f-r,e+r,b['mask'])
            if rect[0]<rect[1] and rect[2]<rect[3]:rs.append(rect)
    return inner,outer


def exact_objective(points,masses,boxes,S,bits):
    inner,outer=exact_rectangles(points,masses,boxes,bits)
    keep=[i for i in range(len(points)) if not S>>i&1]
    pp=[points[i] for i in keep];mm=[masses[i] for i in keep]
    rb=closures(pp);ri,ro=exact_rectangles(pp,mm,rb,bits)
    lo=difference_area(inner,ro);hi=difference_area(outer,ri)
    assert 0<=lo<=hi
    D2=1<<(2*bits)
    return {'lower':str(F(lo,D2)),'upper':str(F(hi,D2)),'bits':bits,
            'full_inner_area':str(F(difference_area(inner,[]),D2)),
            'full_outer_area':str(F(difference_area(outer,[]),D2)),
            'remaining_inner_area':str(F(difference_area(ri,[]),D2)),
            'remaining_outer_area':str(F(difference_area(ro,[]),D2)),
            'full_active_boxes':len(outer),'remaining_active_boxes':len(ro),
            'scope':'Exact rational enclosure of complete true H(full)\\H(remaining); S mass exactly1.'}


def float_area(points,masses,boxes):
    rects=[]
    for b in boxes:
        mass=sum(masses[i] for i in b['ids'])
        if mass<=float(b['span2']):continue
        r=math.sqrt(mass)/2;a,c,e,f=map(float,b['bbox'])
        rects.append((c-r,a+r,f-r,e+r))
    return difference_area(rects,[])


def proposal_objective(points,masses,boxes,keep,residual_boxes):
    return float_area(points,masses,boxes)-float_area([points[i] for i in keep],[masses[i] for i in keep],residual_boxes)


def atomic_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--bits',type=int,default=72)
    ap.add_argument('--trust-radius',default='1/16');ap.add_argument('--finite-diff-step',default='1/262144')
    ap.add_argument('--rounding-bits',type=int,default=28);ap.add_argument('--line-steps',type=int,default=8)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert args.bits>=16 and args.rounding_bits>=16 and 1<=args.line_steps<=16
    import numpy as np
    import scipy
    from scipy.optimize import linprog
    started=time.perf_counter();spec=json.loads(args.input.read_text())
    points=[tuple(map(F,p)) for p in spec['points']];seed=list(map(F,spec['initial_masses']));n=len(points)
    boxes=closures(points);rows=constraints(spec,boxes);check=exact_check(seed,spec,rows,boxes)
    assert check['feasible'],check
    S=spec['S_mask'];keep=[i for i in range(n) if not S>>i&1];residual_boxes=closures([points[i] for i in keep])
    original=exact_objective(points,seed,boxes,S,args.bits)
    payload={'schema':'r064-one-fixed-position-mass-LP-v1','status':'PROPOSAL_STARTED',
             'input_sha256':sha(args.input),'source_sha256':sha(Path(__file__)),
             'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
             'parameters':vars(args)|{'input':str(args.input),'out':str(args.out)},
             'scope':'One fixed21-source input. Floating finite differences/LP propose only; no certified chamber gradient. Exact rational feasibility + complete halo interval decide improvement.',
             'closed_source_masks':len(boxes),'seed_check':check,'seed_objective':original,'candidates':[]}
    atomic_json(args.out,payload)
    # Tangent coordinate proposals preserve sum_S=1. At symmetric walls these
    # are finite differences, not an asserted ordinary derivative.
    epsilon=float(F(args.finite_diff_step));s0=list(map(float,seed));anchor=spec['S_indices'][-1]
    gradient=np.zeros(n)
    for i in range(n):
        if i==anchor:continue
        plus=s0.copy();minus=s0.copy();plus[i]+=epsilon;minus[i]-=epsilon
        if S>>i&1:plus[anchor]-=epsilon;minus[anchor]+=epsilon
        assert min(plus)>0 and min(minus)>0
        gradient[i]=(proposal_objective(points,plus,boxes,keep,residual_boxes)-proposal_objective(points,minus,boxes,keep,residual_boxes))/(2*epsilon)
    A=[];b=[]
    for mask,rhs,kind in rows:
        sign=-1 if kind=='negative-connection-lower' else 1
        A.append([sign*int(mask>>i&1) for i in range(n)]);b.append(float(rhs))
    floor=F(spec['mass_constraints']['positive_floor_for_initial_chamber']);trust=F(args.trust_radius)
    assert trust>0
    bounds=[(float(max(floor,m-trust)),float(m+trust)) for m in seed]
    eq=[[int(S>>i&1) for i in range(n)]]
    lp=linprog(-gradient,A_ub=np.asarray(A),b_ub=np.asarray(b),A_eq=np.asarray(eq),b_eq=[1.0],bounds=bounds,method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    payload['gradient_proposal']=gradient.tolist();payload['LP']={'success':bool(lp.success),'status':lp.status,'message':lp.message,'predicted_gain':None if lp.fun is None else float(-lp.fun-gradient@np.asarray(s0))}
    if not lp.success:
        payload['status']='LP_FAILED_NO_MATHEMATICAL_CLAIM';atomic_json(args.out,payload);return 2
    D=1<<args.rounding_bits
    vertex=[F(round(float(x)*D),D) for x in lp.x]
    vertex[anchor]=1-sum(vertex[i] for i in spec['S_indices'] if i!=anchor)
    direction=[v-m for v,m in zip(vertex,seed)]
    payload['rational_proposal_vertex']=list(map(str,vertex));payload['rational_direction']=list(map(str,direction))
    for k in range(args.line_steps):
        step=F(1,1<<k);candidate=[m+step*v for m,v in zip(seed,direction)]
        c=exact_check(candidate,spec,rows,boxes);record={'step':str(step),'masses':list(map(str,candidate)),'exact_check':c}
        if c['feasible']:
            record['objective']=exact_objective(points,candidate,boxes,S,args.bits)
            record['certified_strict_improvement']=F(record['objective']['lower'])>F(original['upper'])
            record['certified_gain_lower']=str(F(record['objective']['lower'])-F(original['upper']))
        payload['candidates'].append(record);atomic_json(args.out,payload)
    winners=[r for r in payload['candidates'] if r.get('certified_strict_improvement')]
    payload['status']='CERTIFIED_INPUT_SPECIFIC_IMPROVEMENT' if winners else 'NO_CERTIFIED_IMPROVEMENT_ON_THIS_LINE'
    if winners:payload['best_candidate']=max(winners,key=lambda r:F(r['objective']['lower']))
    payload['wall_seconds']=time.perf_counter()-started;atomic_json(args.out,payload)
    print(json.dumps({'status':payload['status'],'candidate_count':len(payload['candidates']),'seconds':payload['wall_seconds']}))
    return 0


if __name__=='__main__':sys.exit(main())
