#!/usr/bin/env python3
"""Certify two consecutive hexagon poses, choosing by the complementary cap."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import audit_optimized_eh_chord as geo
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


def certify(program,right,max_depth,replay=None):
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
    for turn,complementary_cap in ((0,'K-C3'),(1,'I-E3')):
        for label,normal,offset in inner:
            if label in hex_labels or label==complementary_cap:
                continue
            key=f'{turn}:{label}'
            cache={}
            used=set()
            cells=[]
            local=[]

            def end(angle,mu):
                location=f'{angle}_{mu}'
                if location not in cache:
                    objective,threshold=root_clearance_objective(program,angle-60*turn,normal,offset,mu)
                    weights=program.optimize(objective) if replay is None else replay['edges'][key]['weights'][location]
                    cache[location]=(program.validate(objective,weights),threshold,weights)
                return location,cache[location]

            def margin(a,b,mu):
                ka,(ua,threshold,_)=end(a,mu)
                kb,(ub,_,_)=end(b,mu)
                cosine=(exact_fraction(b-a)*geo.PI/360).cos()
                assert cosine>0 and max(ua,ub)>0
                return threshold-max(ua,ub)/cosine,ka,kb

            def record(a,b,mu,clearance,ka,kb):
                assert clearance>0
                cells.append([str(a),str(b),str(mu)])
                used.update((ka,kb))
                local.append(clearance.lower())

            def split(a,b,depth):
                for mu in multipliers(label):
                    clearance,ka,kb=margin(a,b,mu)
                    if clearance>0:
                        record(a,b,mu,clearance,ka,kb)
                        return
                if depth==max_depth:
                    raise ArithmeticError(f'uncertified {key} on {a}--{b}')
                center=(a+b)/2
                split(a,center,depth+1)
                split(center,b,depth+1)

            if replay is None:
                split(left,right,0)
            else:
                previous=left
                for aa,bb,mm in replay['edges'][key]['cells']:
                    a,b,mu=map(Fraction,(aa,bb,mm))
                    assert left<=a==previous<b<=right
                    record(a,b,mu,*margin(a,b,mu))
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
    return {'summary':result,'edges':results,'right_weights':right_weights,'adaptive_cuts':program.adaptive_cuts}


def run(args):
    sigma=Fraction(args.sigma)
    constructor,model_hash=specialize(sigma)
    trial=geo.trial(args.delta,tight=True)
    signs=voronoi_signs(arc_data())
    program=constructor(args.delta)
    right=Fraction(args.right)
    replay=None
    if args.replay:
        replay=json.loads(args.replay.read_text())
        source=replay['summary']
        assert source['delta']==args.delta and Fraction(source['sigma_degrees'])==sigma
        assert Fraction(source['root_angle_degrees'][1])==right
        assert source['specialized_model_sha256']==model_hash
        for index,g in replay['adaptive_cuts']:
            program.cut(index,g,record=True)
    output=certify(program,right,args.max_depth,replay)
    summary=output['summary']
    summary['sigma_degrees']=str(sigma)
    summary['specialized_model_sha256']=model_hash
    summary['Sprague_endpoint_signs']=signs
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
    parser.add_argument('--old-area',type=Path)
    parser.add_argument('--replay',type=Path)
    run(parser.parse_args())
