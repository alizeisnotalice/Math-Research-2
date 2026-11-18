"""Certify a third EH face meeting a deeper point of the true C3 arc."""
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import certify_eh_dual_hinge as prev
from certify_eh_general_relocation import exact_fraction as arbq

eh, g = prev.eh, prev.g
DEPTH = '0.000557'
PIVOT = F('0.5')
RIGHT = F('36.01')
F_LIMIT = F('-0.00000000001')


def unit_outer(a,b):
    direction = g.minus(b,a)
    return g.scale(-1/g.norm(direction),g.rot(direction))


def geometry(old):
    normals,offsets,old_gain,u,w=eh.clip_geometry(old['parameters'])
    reduced=[h-arbq(F(s)) for h,s in zip(offsets,prev.RELAX)]
    x=prev.point_on_faces(normals,reduced)
    t=g.setup()[7]
    y=g.plus(x,g.scale(arbq(PIVOT),g.minus(u,x)))
    deeper=eh.clip_geometry([DEPTH,*old['parameters'][1:]])[3]
    new_face=unit_outer(deeper,y)
    normals=[*normals,new_face]; reduced=[*reduced,g.dot(new_face,deeper)]
    points=[deeper,y,x,t];edges=[g.minus(b,a) for a,b in zip(points,points[1:])]
    assert all(v[1]<0 for v in edges) and all(g.cross(a,b)>0 for a,b in zip(edges,edges[1:]))
    _,_,c3,b3,*_=g.setup()
    assert g.cross(g.rot(g.minus(deeper,c3)),edges[0])>0
    assert g.cross(edges[-1],g.rot(g.minus(t,b3)))>0
    assert g.cross(edges[0],g.minus(u,deeper))<0
    assert all(g.cross(g.minus(t,deeper),g.minus(a,deeper))<0 for a in (y,x))
    # The new top facet and the old upper facet meet at y. On the lower side
    # of that height, the new top facet is redundant by the turn sign.
    assert g.cross(edges[0],edges[1])>0
    trial=g.trial(DEPTH,tight=True)
    turn_loss=sum(g.cross(g.minus(b,deeper),g.minus(c,deeper)) for b,c in zip(points[1:-1],points[2:]))/2
    gain=g.arb(trial['gain_area'])-turn_loss-old_gain
    assert gain>0
    return normals,reduced,points,u,gain


def branch_model(model,order,normals,offsets,height):
    deep='2' in order
    p=model(DEPTH if deep else '0.000555')
    basis=np.eye(p.dimension+1,dtype=object)
    for i,letter in enumerate(order):
        j=int(letter);point=basis[2*i:2*i+2]
        p.linear(-np.asarray(normals[j],dtype=object)@point,-offsets[j])
        if j==2:
            p.linear(-point[1],-height)
    return p


def replay_cuts(program,cuts):
    for index,direction in cuts:
        program.cut(index,direction,record=True)


def trim(payload,weights,initial):
    # Only the elementary supporting inequalities carrying nonzero weight are needed.
    active=sorted({i-initial for row in weights for i,_ in row if i>=initial})
    kept={initial+i:initial+j for j,i in enumerate(active)}
    payload['adaptive_cuts']=[payload['adaptive_cuts'][i] for i in active]
    for row in weights:
        row[:]=[[kept.get(i,i),w] for i,w in row]


def run(archive,ah_path,seed,emit):
    path=archive/'eh-two-line-cover-certificate.json';old=json.loads(path.read_bytes());ah=json.loads(ah_path.read_bytes())
    assert old['parameters']==['0.000555','0.0005522','0.0075']
    if 'source_hashes' in ah:
        assert ah['source_hashes'][path.name]==sha256(path.read_bytes()).hexdigest()
    else:
        assert ah['R167_area_sha256']==sha256((archive/'ah-three-line-cut-certificate.json').read_bytes()).hexdigest()
    model,model_hash=eh.specialize(F('1.5497')); assert model_hash==old['specialized_model_sha256']
    normals,offsets,points,u,gain=geometry(old);deeper,y,x,t=points
    fresh=seed is None
    if fresh: seed={'new_branches':{}}
    f_bounds={}
    for order in ('00','01','10','11','20','21','02','12','22'):
        p=branch_model(model,order,normals,offsets,y[1]);initial=len(p.rows);obj=eh.phi_objective(p,RIGHT,1)
        if '2' not in order:
            s=old['ordered_branches'][order];replay_cuts(p,s['adaptive_cuts']);weights=s['dual']
        elif fresh:
            for s in old['ordered_branches'].values():replay_cuts(p,s['adaptive_cuts'])
            weights=p.optimize(obj);s={'adaptive_cuts':list(p.adaptive_cuts),'dual':deepcopy(weights)}
            trim(s,[s['dual']],initial);seed['new_branches'][order]=s
        else:
            s=seed['new_branches'][order];replay_cuts(p,s['adaptive_cuts']);weights=s['dual']
        f=p.validate(obj,weights)-g.arb(3)/2
        assert f<arbq(F_LIMIT),(order,str(f));f_bounds[order]=str(f.upper())
    p=model(DEPTH);eh.add_joint_query(p);p.linear(eh.phi_objective(p,RIGHT,1),g.arb(3)/2+arbq(F_LIMIT));initial=len(p.rows)
    if fresh:
        stencil=deepcopy(old['enclosing_pose_certificate']);replay_cuts(p,stencil['adaptive_cuts']);validate=p.validate
        def install(obj,weight):
            if weight is not stencil['right_weights']:
                weight[:]=p.optimize(obj)
            return validate(obj,weight)
        p.validate=install
        pose=eh.certify(p,RIGHT,10,stencil,eh_relocation=True,joint_eh=True,revisit_ah=True)
        weights=[pose['right_weights']]+[row for group in pose['edges'].values() for row in group['weights'].values()]
        trim(pose,weights,initial);pose.pop('summary');seed['pose']=pose
        emit.write_text(json.dumps(seed,separators=(',',':'))+'\n')
        return
    stencil=seed['pose'];replay_cuts(p,stencil['adaptive_cuts']);pose=eh.certify(p,RIGHT,10,stencil,eh_relocation=True,joint_eh=True,revisit_ah=True)['summary']
    old_n,new_n=unit_outer(u,t),unit_outer(deeper,t)
    # Normal-change allowance over the entire hexagon, using the common endpoint T.
    cost=2*g.norm(g.minus(new_n,old_n))/g.arb(3).sqrt()
    ah_margin=min(g.arb(row['placement_margin_lower']).lower() for row in ah['bins']);assert ah_margin>cost
    area_low,area_high=map(F,ah['area_exact_interval']);lo,hi=map(eh.gain_endpoint,(gain.lower(),gain.upper()))
    bounds=area_low-hi,area_high-lo
    return {'depth':DEPTH,'hinge_fraction_exact':str(PIVOT),'F_limit_exact':str(F_LIMIT),'model_sha256':model_hash,
            'AH_input_sha256':sha256(ah_path.read_bytes()).hexdigest(),'R126_input_sha256':sha256(path.read_bytes()).hexdigest(),
            'ordered_F_upper':f_bounds,'EH_placement':pose,'boundary_vertices':[list(map(str,v)) for v in points],
            'original_AH_support_margin':str(ah_margin),'AH_chord_change_cost_upper':str(cost.upper()),'AH_relocated_support_margin_lower':str((ah_margin-cost).lower()),
            'minimum_AH_clearance_whole_EH_rectangle':prev.ah_clearance(ah,[(a,b) for a in (deeper[0],t[0]) for b in (deeper[1],t[1])]),
            'sprague_endpoint_signs':eh.voronoi_signs(eh.arc_data()),'EH_gain_over_R126_exact_interval':[str(lo),str(hi)],
            'gain_over_R126_arb':str(gain),'area_exact_interval':list(map(str,bounds)),'area_decimal_interval':prev.format_area(*bounds)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('archive',type=Path);parser.add_argument('AH',type=Path);parser.add_argument('weights',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    if args.generate:run(args.archive,args.AH,None,args.weights)
    result=run(args.archive,args.AH,json.loads(args.weights.read_bytes()),None);result['weights_sha256']=sha256(args.weights.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
