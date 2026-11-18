#!/usr/bin/env python3
"""Certify two consecutive hexagon poses, choosing by the complementary cap."""
import argparse
from fractions import Fraction
import json
from math import ceil
from pathlib import Path

import numpy as np

import audit_optimized_eh_chord as geo
from certify_eh_general_relocation import numeric, rotate_angle
from certify_eh_joint_parameters import gain_endpoint, specialize
from certify_eh_relocation_root_condition import (
    exact_fraction, phi_objective, root_clearance_objective,
)
from certify_eh_sprague_relocation import arc_data, hybrid_polygon, voronoi_signs


def multipliers(label):
    if label=='U1-Te':
        choices=('0','-.128','-.125','-.13','-.12','-.135','-.14','-.2','-.1','-.15','-.25','-.3','-.5','1','-1','2','-2')
    else:
        choices=('0','-1','-.5','1','1.5','-1.5','2','-2','3','-3','5')
    return tuple(map(Fraction,choices))


def add_joint_query(program):
    old=program.dimension
    program.dimension+=2
    eye=np.eye(program.dimension+1,dtype=object)

    def extend(p):
        return np.column_stack((p[:,:old],np.zeros((2,2),dtype=object),p[:,-1]))

    program.witnesses=[extend(p) for p in program.witnesses]
    program.support=[extend(p) for p in program.support]
    program.disks=[extend(p) for p in program.disks]
    program.rows=[(np.concatenate((a,np.zeros(2,dtype=object))),b) for a,b in program.rows]
    program.numeric_rows=[numeric(a) for a,_ in program.rows]
    extra=eye[old:old+2]
    program.support.append(extra)
    n,m,*_=geo.setup()
    axis=n[5]
    reflection=2*np.outer(axis,axis)-np.eye(2,dtype=object)
    for u in n+[m[3],m[5],reflection@m[3],reflection@m[5]]:
        program.linear(np.asarray(u,dtype=object)@extra,geo.HALF)
    first=len(program.disks)
    program.disks.extend(extra-y for y in program.witnesses)
    program.disk_num=np.stack([numeric(d) for d in program.disks])
    for index in range(first,len(program.disks)):
        for direction in ((1,0),(-1,0),(0,1),(0,-1)):
            program.cut(index,direction)


def joint_objective(program,angle,normal,image,offset,multiplier,ratio):
    mu,a=map(exact_fraction,(multiplier,ratio))
    average=geo.plus(geo.scale(a,normal),geo.scale(1-a,image))
    root3=geo.arb(3).sqrt()
    d0,d1=average[0]-average[1]/root3,2*average[1]/root3
    directions=[normal,image,*program.normals[:3]]
    weights=[a,1-a,mu-d0,-mu-d1,mu]
    threshold=offset-(d0+d1-mu)/2
    for i,w in enumerate(weights):
        if w<0:
            weights[i]=-w
            directions[i]=geo.neg(directions[i])
            threshold-=w
        else:
            assert w>0 or w.is_zero()
    rotation=exact_fraction(angle)*geo.PI/180-4*geo.PI/3
    directions=[rotate_angle(v,rotation) for v in directions]
    return program.objective(directions,weights),threshold


def certify(program,right,max_depth,replay=None,eh_relocation=False,joint_eh=False,revisit_ah=False):
    left=Fraction(30)
    objective=phi_objective(program,right,1)
    right_weights=program.optimize(objective) if replay is None else replay['right_weights']
    endpoint_upper=program.validate(objective,right_weights)-geo.arb(3)/2
    assert endpoint_upper<0,str(endpoint_upper)
    inner,_=hybrid_polygon(program.delta,4)
    hex_labels={'B1-Q','P-F1','F1-H','E3-D1','D1-G','C3-B1'}
    assert len([edge for edge in inner if edge[0] in hex_labels])==6
    assert {'I-E3','K-C3'}<={label for label,_,_ in inner}
    results={}
    margins=[]
    reflective_cells=0
    joint_cells=0
    axis=program.normals[5]

    def reflect(u):
        return geo.minus(geo.scale(2*geo.dot(axis,u),axis),u)

    for turn,complementary_cap in ((0,'K-C3'),(1,'I-E3')):
        for label,normal,offset in inner:
            if label in hex_labels or label==complementary_cap:
                continue
            if revisit_ah and label not in ('I-E3','K-C3','U1-Te'):
                continue
            key=f'{turn}:{label}'
            cache={}
            used=set()
            cells=[]
            local=[]

            def end(angle,mu,ratio):
                if ratio==1:
                    location=f'{angle}_{mu}'
                elif ratio==0:
                    location=f'{angle}_{mu}_Q2'
                else:
                    location=f'{angle}_{mu}_Q2mean={ratio}'
                if location not in cache:
                    if ratio in (0,1):
                        measured=normal if ratio==1 else reflect(normal)
                        objective,threshold=root_clearance_objective(program,angle-60*turn,measured,offset,mu)
                    else:
                        objective,threshold=joint_objective(program,angle-60*turn,normal,reflect(normal),offset,mu,ratio)
                    weights=program.optimize(objective) if replay is None else replay['edges'][key]['weights'][location]
                    cache[location]=(program.validate(objective,weights),threshold,weights)
                return location,cache[location]

            def margin(a,b,mu,ratio):
                ka,(ua,threshold,_)=end(a,mu,ratio)
                kb,(ub,_,_)=end(b,mu,ratio)
                cosine=(exact_fraction(b-a)*geo.PI/360).cos()
                assert cosine>0 and max(ua,ub)>0
                return threshold-max(ua,ub)/cosine,ka,kb

            def record(a,b,mu,ratio,clearance,ka,kb):
                nonlocal reflective_cells,joint_cells
                assert clearance>0
                row=[str(a),str(b),str(mu)]
                if ratio==0:
                    assert eh_relocation and label=='U1-Te'
                    row.append('Q2')
                    reflective_cells+=1
                elif ratio!=1:
                    assert joint_eh and eh_relocation and label=='U1-Te' and 0<ratio<1
                    row.append(f'Q2mean={ratio}')
                    joint_cells+=1
                cells.append(row)
                used.update((ka,kb))
                local.append(clearance.lower())

            def split(a,b,depth):
                ratios=[Fraction(1)]
                if eh_relocation and label=='U1-Te':
                    ratios.append(Fraction(0))
                    if joint_eh:
                        ratios.extend((Fraction(1,4),Fraction(1,2),Fraction(3,4)))
                for ratio in ratios:
                    base=multipliers(label)
                    if label=='U1-Te' and joint_eh:
                        if ratio==1:
                            base=tuple(map(Fraction,('0','-.08','-.07','-.065','-.1','-.13')))+base
                        elif ratio==0:
                            base=tuple(map(Fraction,('0','.01','.015')))+base
                        else:
                            base=tuple(map(Fraction,('-.015','-.03','-.05','-.025','0','-.04','-.055','-.07')))
                    for mu in base:
                        clearance,ka,kb=margin(a,b,mu,ratio)
                        if clearance>0:
                            record(a,b,mu,ratio,clearance,ka,kb)
                            return
                if depth==max_depth:
                    raise ArithmeticError(f'uncertified {key} on {a}--{b}')
                center=(a+b)/2
                split(a,center,depth+1)
                split(center,b,depth+1)

            if replay is None:
                presets=[]
                if joint_eh and label=='U1-Te' and turn==0 and right>31:
                    presets=[('30','30.1','0','.01'),('30.1','30.25','1/4','-.015'),
                             ('30.25','30.4','3/4','-.04'),('30.4','31','1','-.05')]
                for preset in presets:
                    start,finish,ratio,mu=map(Fraction,preset)
                    value,ka,kb=margin(start,finish,mu,ratio)
                    if value>0:
                        record(start,finish,mu,ratio,value,ka,kb)
                    else:
                        split(start,finish,0)
                start=Fraction(presets[-1][1]) if presets else left
                if revisit_ah:
                    for stop in range(61,ceil(2*right)):
                        endpoint=Fraction(stop,2)
                        if start<endpoint<right:
                            split(start,endpoint,0)
                            start=endpoint
                split(start,right,0)
            else:
                previous=left
                for row in replay['edges'][key]['cells']:
                    assert len(row)==3 or (len(row)==4 and (row[3]=='Q2' or row[3].startswith('Q2mean=')))
                    a,b,mu=map(Fraction,row[:3])
                    ratio=Fraction(1)
                    if len(row)==4:
                        ratio=Fraction(0) if row[3]=='Q2' else Fraction(row[3].split('=',1)[1])
                    assert left<=a==previous<b<=right
                    record(a,b,mu,ratio,*margin(a,b,mu,ratio))
                    previous=b
                assert previous==right
            results[key]={'cells':cells,'weights':{place:cache[place][2] for place in sorted(used)}}
            lowest=min(local)
            margins.append((lowest,key))
            print(key,len(cells),str(lowest),flush=True)

    worst,key=min(margins,key=lambda row:row[0])
    result={'delta':program.delta,'root_angle_degrees':[str(left),str(right)],
            'F_right_upper':str(endpoint_upper),'test_edge_count':len(results),'active_edge':key,
            'minimum_slack_lower_bound':str(worst),'adaptive_cut_rows':len(program.adaptive_cuts),
            'angle_edge_cells':sum(len(row['cells']) for row in results.values()),
            'validated_duals':1+sum(len(row['weights']) for row in results.values())}
    if eh_relocation:
        result['EH_relocation_cells']=reflective_cells
    if joint_eh:
        result['EH_combined_cells']=joint_cells
    if revisit_ah:
        result['AH_relocation_after_root']=True
    return {'summary':result,'edges':results,'right_weights':right_weights,'adaptive_cuts':program.adaptive_cuts}


def run(args):
    sigma=Fraction(args.sigma)
    constructor,model_hash=specialize(sigma)
    trial=geo.trial(args.delta,tight=True)
    signs=voronoi_signs(arc_data())
    program=constructor(args.delta)
    if args.revisit_ah:
        assert args.eh_relocation
    if args.joint_eh:
        assert args.eh_relocation
        add_joint_query(program)
    if args.generation_tolerance is not None:
        assert args.replay is None and float(args.generation_tolerance)>0
        optimize=program.optimize
        program.optimize=lambda objective:optimize(objective,tolerance=float(args.generation_tolerance))
    right=Fraction(args.right)
    replay=None
    junction=None
    if args.eh_relocation:
        _,_,c3,b3,_,_,_,_,j=geo.setup()
        first,second=geo.minus(j,c3),geo.minus(j,b3)
        turn,advance=geo.cross(first,second),geo.dot(first,second)
        assert turn>0 and advance>0
        junction=[str(turn.lower()),str(advance.lower())]
    if args.replay:
        replay=json.loads(args.replay.read_text())
        source=replay['summary']
        assert source['delta']==args.delta and Fraction(source['sigma_degrees'])==sigma
        assert Fraction(source['root_angle_degrees'][1])==right
        assert source['specialized_model_sha256']==model_hash
        assert ('EH_relocation_cells' in source)==args.eh_relocation
        assert ('EH_combined_cells' in source)==args.joint_eh
        assert ('AH_relocation_after_root' in source)==args.revisit_ah
        for index,g in replay['adaptive_cuts']:
            program.cut(index,g,record=True)
    output=certify(program,right,args.max_depth,replay,args.eh_relocation,args.joint_eh,args.revisit_ah)
    summary=output['summary']
    summary['sigma_degrees']=str(sigma)
    summary['specialized_model_sha256']=model_hash
    summary['Sprague_endpoint_signs']=signs
    if junction is not None:
        summary['restored_ES_junction_cross_and_dot_lower']=junction
    gain=geo.arb(trial['gain_area'])
    lo,hi=(gain_endpoint(point) for point in (gain.lower(),gain.upper()))
    summary['cut_gain_area_interval_exact']=[str(lo),str(hi)]
    if args.old_area:
        old=json.loads(args.old_area.read_text())
        assert Fraction(old['sigma_degrees'])==sigma
        a,b=map(Fraction,old['area_interval'])
        summary['candidate_area_interval_exact']=[str(a-hi),str(b-lo)]
    args.output.write_text(json.dumps(output,separators=(',',':'))+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sigma')
    parser.add_argument('delta')
    parser.add_argument('right')
    parser.add_argument('output',type=Path)
    parser.add_argument('--max-depth',type=int,default=8)
    parser.add_argument('--eh-relocation',action='store_true')
    parser.add_argument('--joint-eh',action='store_true')
    parser.add_argument('--revisit-ah',action='store_true')
    parser.add_argument('--generation-tolerance')
    parser.add_argument('--old-area',type=Path)
    parser.add_argument('--replay',type=Path)
    run(parser.parse_args())
