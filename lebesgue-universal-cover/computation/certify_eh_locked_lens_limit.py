"""Bound arbitrary convex refinements using two uniquely placed constant-width completions."""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import certify_eh_upper_chain as chain
from certify_eh_no_interior_roots import verify as verify_rigidity
from certify_eh_polyline_obstruction import witness_points

g=chain.g


def mix(source, index, coefficient):
    anchor=source['hinge_interval_witnesses'][9];other=source['new_witnesses'][index]
    assert anchor['support_angle_degrees']==other['support_angle_degrees']=='36.01' and anchor['sigma']==other['sigma']=='1.5497'
    with localcontext() as ctx:
        ctx.prec=100
        k=Decimal(coefficient)
        a,b=anchor['coordinates'],other['coordinates']
        coordinates=[format((1-k)*Decimal(x)+k*Decimal(y),'f') for x,y in zip(a,b,strict=True)]
    for x,y,z in zip(a,b,coordinates,strict=True): assert F(z)==(1-F(coefficient))*F(x)+F(coefficient)*F(y)
    return {'sigma':'1.5497','delta':'.001','support_angle_degrees':'36.01','coordinates':coordinates}


def verify(archive,parameter_path):
    source_path=archive/'eh-polyline-obstruction-input.json';source=json.loads(source_path.read_bytes());params=json.loads(parameter_path.read_bytes());old=json.loads((archive/'eh-two-line-cover-certificate.json').read_bytes())
    chain.b.eh.specialize(F('1.5497'));front=chain.front(old,params);v,t=front['v'],front['t'];core_normal=chain.b.unit_outer(v,t)
    settings=[(1,'.000006019836',False,('b','contact_0'),1), (5,'-.0003784',True,('b_Q','contact_0_op'),-1)]
    fixed=[];certs=[]
    for index,weight,reflected,pair,side in settings:
        case=mix(source,index,weight);receipt=verify_rigidity(case)
        assert receipt['num_witness_points']==18 and receipt['only_d6_poses_in_P']==['identity','Q2'] and len(receipt['excluded_d6_poses'])==10
        points,Q=witness_points(case);outside=points['z'] if reflected else Q(points['contact_0'])
        violation=g.dot(front['normals'][1],outside)-front['offsets'][1];assert violation>0
        a,bb=(points[k] for k in pair);dist=g.norm(g.minus(a,bb));assert 0<dist<2
        extra=g.circle_pair(a,bb,side)
        remaining=min((1-g.dot(g.minus(extra,pt),g.minus(extra,pt))).lower() for k,pt in points.items() if k not in pair);assert remaining>0
        forced=Q(extra) if reflected else extra;chord_distance=g.dot(core_normal,g.minus(forced,t));assert chord_distance>0
        fixed.append(forced);certs.append({'source_interpolation_index':index,'coefficient':weight,'new_coordinates':case['coordinates'],'allowed_pose':'Q2' if reflected else 'identity',
            'excluded_other_pose_margin_lower':str(violation.lower()),'two_active_unit_circles':list(pair),'intersection_sign':side,'other_sixteen_squared_slack_lower':str(remaining),
            'forced_point_outside_core_distance_lower':str(chord_distance.lower()),'full_rigidity_replay':{k:receipt[k] for k in ('squared_slack_other_pairs','uniform_F_prime_lower_first_quarter_degree','uniform_F_lower_interior','interior_angle_cell_count')}})
    tri=lambda a,b,c:g.cross(g.minus(b,a),g.minus(c,a))/2
    a,bb=fixed;boundary=[v,a,bb,t];turns=[tri(boundary[i],boundary[(i+1)%4],boundary[(i+2)%4]) for i in range(4)];assert min(turns)>0
    forced_area=tri(v,a,bb)+tri(v,bb,t);assert forced_area>0
    core_loss=g.arb(g.trial(params[-1][0],tight=True)['gain_area'])-chain.b.eh.clip_geometry(old['parameters'])[2]-front['gain']
    headroom=core_loss-forced_area;assert 0<headroom<g.arb('8.645331e-12');return {'witness_source_sha256':sha256(source_path.read_bytes()).hexdigest(),'profile_sha256':sha256(parameter_path.read_bytes()).hexdigest(),
        'completions':certs,'four_positive_oriented_half_determinants':list(map(str,turns)),'mandatory_quadrilateral_area':str(forced_area),
        'EH_chain_minus_core_area':str(core_loss),'refinement_gain_upper_exact':str(chain.b.eh.gain_endpoint(headroom.upper())),'refinement_gain_upper':str(headroom.upper())}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('parameters',type=Path);p.add_argument('output',type=Path);s=p.parse_args();s.output.write_text(json.dumps(verify(s.archive,s.parameters),indent=2)+'\n')
