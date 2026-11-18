"""Certify a positive-angular interior box by trading one weight for y3 motion."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb,ctx
from certify_three_body_twist_box import (
    base_topology, chain, contacts, jets, other_root, pi, q, shifted_bound,
    short_area, solve3, solve_interval, subintervals, transpose, transverse_bounds,
)
from certify_three_body_symmetric_slice import interval

A=('-.017','-.0076')
PHI_BEGIN=F(147,32768)
PHI_END=F(1,50)
OTHER=(('-.0008','.0008'),('.0185','.0274'),('-.0045','.001'))
COEFFICIENT=F(9,100)
BASELINE_CUTS=('-.0147','-.0098')
M0_UPPER=q('.8347809457752037')


def margins(a,phi,upper,actual,full_y,check_gain,baseline):
    y,ey=solve_interval(a,phi,upper,full_y)
    _,_,coef,T1=jets(a,y,phi,upper)
    result={'implicit_root':ey.lower(), 'curvature_cushion':(coef-q(COEFFICIENT)).lower(),
            'upper_angular_slope':(1-abs(T1)).lower()}
    if not check_gain:
        return result
    _,data=contacts(a,y,phi,upper)
    _,_,_,_,_,wplus,wminus,np,nm,nT,t5=data
    parts,angles=chain(data,phi,upper)
    gp,dp,gm,dm,T=angles
    result['angle_branch_and_layout']=min(x.lower() for x in [
        np[0],nm[0],other_root(wplus,np)[0],other_root(wminus,nm)[0],nT[1],-wplus[0],-wminus[0],
        gp-pi/10-phi,pi/6-gp,gm+pi/6,-pi/10+phi-gm,dp-pi/6,3*pi/10+phi-dp,
        dm+3*pi/10-phi,-pi/6-dm,T-pi*q('.49'),pi*q('.51')-T])
    weights,_,B=transverse_bounds(data,parts,upper)
    move=solve3(transpose(B),[1,0,0])
    assert move[1]>8
    correction=max(arb(0),((weights[1]-q('.99999')).upper()/move[1].lower()).upper())
    legal=[w-correction*v for w,v in zip(weights,move)]
    result['weight_after_y_adjustment']=min((1-abs(v)).lower() for v in legal)
    result['fan_connector_determinant']=shifted_bound(data,parts,actual,upper)[0]

    dy=actual[0]-y
    db=actual[1]-t5[0]
    nonlinear=abs(legal[1]-legal[0]-legal[2]).abs_upper()*abs(dy*db).abs_upper()/2
    gain=baseline+q(COEFFICIENT)*phi.lower()**2+correction*dy.lower()-nonlinear
    result['gain_remaining_after_transverse_adjustment']=gain.lower()
    return result


def main():
    parser=ArgumentParser()
    parser.add_argument('--a-count',type=int,default=470)
    parser.add_argument('--early-count',type=int,default=33)
    parser.add_argument('--angle-count',type=int,default=112)
    parser.add_argument('--a-start',type=int,default=0)
    parser.add_argument('--a-stop',type=int)
    args=parser.parse_args()
    assert args.a_count>0 and args.early_count>0 and args.angle_count>0
    a_end=args.a_stop if args.a_stop is not None else args.a_count
    assert 0<=args.a_start<a_end<=args.a_count
    early=list(subintervals(0,PHI_BEGIN,args.early_count))
    later=list(subintervals(PHI_BEGIN,PHI_END,args.angle_count))
    actual=tuple(interval(*pair) for pair in OTHER)
    full_y=interval('-.001','.001')
    tracked={}
    where={}
    ranges=list(subintervals(*A,args.a_count))
    cut1,cut2=map(q,BASELINE_CUTS)
    gaps=[]
    for cut in (cut1,cut2):
        details=contacts(cut,arb(0),arb(0),True)[1]
        value=short_area(details,arb(0),True)[0]
        assert value>M0_UPPER
        gaps.append((value-M0_UPPER).lower())
    for i,a in enumerate(ranges[args.a_start:a_end],args.a_start):
        assert base_topology(a)>0
        base_b=contacts(a,arb(0),arb(0),True)[1][-1][0]
        assert q('.0185')<base_b<q('.0274')
        baseline=gaps[0] if a<cut1 else gaps[1] if a>cut2 else arb(0)
        for j,phi in enumerate(early+later):
            for upper in (True,False):
                values=margins(a,phi,upper,actual,full_y,j>=len(early),baseline)
                for key,val in values.items():
                    assert val>0,(key,i,j,upper,val,a,phi)
                    if key not in tracked or val<tracked[key]:
                        tracked[key]=val
                        where[key]=[i,j,upper]
    print(json.dumps({'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'base_scope':[str(0),str(PHI_END)],'angle_gain_scope':[str(PHI_BEGIN),str(PHI_END)],
        'box_a':A,'box_y_b_z':OTHER, 'a_partition':args.a_count,'a_checked':[args.a_start,a_end],
        'early_count':args.early_count,'angle_count':args.angle_count,
        'angular_quadratic_lower':str(COEFFICIENT), 'baseline_cuts':BASELINE_CUTS, 'baseline_gaps':list(map(str,gaps)),
        'cells':2*(a_end-args.a_start)*(args.early_count+args.angle_count),
        'strict_margins':{k:str(v) for k,v in tracked.items()},'worst':where},indent=2))


if __name__=='__main__':
    main()
