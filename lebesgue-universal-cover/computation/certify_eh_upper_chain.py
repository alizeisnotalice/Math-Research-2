"""Verify an EH chain by one sign inequality per upper facet and four common high-point cases."""
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import certify_eh_third_facet as b
from certify_eh_local_chord_transfer import placement_lower

g=b.g


def front(old,params):
    n,cap,prev_gain,u,w=b.eh.clip_geometry(old['parameters']);t=g.setup()[7]
    relaxed=[h-b.arbq(F(r)) for h,r in zip(cap,b.prev.RELAX)]
    x=b.prev.point_on_faces(n,relaxed)
    normals=[b.unit_outer(u,x),b.unit_outer(x,t)]
    offsets=[g.dot(normals[0],x),g.dot(normals[1],x)]
    _,_,c3,b3,*_=g.setup()
    ys=[]; corners=[]
    assert all(F(d)<F(e) and F(p)<F(q) for (d,p),(e,q) in zip(params,params[1:]))
    for d,p in params:
        assert F('0.000555')<F(d) and 0<F(p)<1
        v=b.eh.clip_geometry([d,*old['parameters'][1:]])[3]
        y=g.plus(x,g.scale(b.arbq(F(p)),g.minus(u,x)))
        direction=b.unit_outer(v,y);normals.append(direction);offsets.append(g.dot(direction,v));ys.append(y[1]);corners.append(v)
        a,bb,cc=g.minus(y,v),g.minus(x,y),g.minus(t,x)
        assert all(edge[1]<0 for edge in (a,bb,cc)) and g.cross(a,bb)>0 and g.cross(bb,cc)>0
        assert g.cross(g.rot(g.minus(v,c3)),a)>0 and g.cross(cc,g.rot(g.minus(t,b3)))>0
        assert g.cross(a,g.minus(u,v))<0 and all(g.cross(g.minus(t,v),g.minus(pt,v))<0 for pt in (x,y))
    v=corners[-1];top=len(normals)-1
    # The scalar x at height y above T increases with the upper depth.
    # Enclosing chord lines pass through T, so a single strict comparison at
    # the shallower endpoint orders the complete portions above T.
    assert all(a[1]<c[1] and g.cross(g.minus(t,c),g.minus(a,c))<0
               for a,c in zip([u]+corners[:-1],corners))
    def slacks(z,skip):return [h-g.dot(a,z) for i,(a,h) in enumerate(zip(normals,offsets)) if i not in skip]
    assert all(c>0 for c in slacks(v,[top])+slacks(t,[1]))
    verts=[(v,[top]),(t,[1])]
    for i,j in combinations(range(len(normals)),2):
        z=b.prev.point_on_faces([normals[i],normals[j]],[offsets[i],offsets[j]])
        if z[1]<t[1] or z[1]>v[1]:continue
        assert t[1]<z[1]<v[1]
        ss=slacks(z,[i,j]);assert all(s>0 or s<0 for s in ss)
        if all(s>0 for s in ss):verts.append((z,[i,j]))
    verts.sort(key=lambda r:r[0][1].mid(),reverse=True)
    edges=[]; active=[]
    for (a,I),(c,J) in zip(verts,verts[1:]):
        common=set(I)&set(J);assert len(common)==1 and c[1]<a[1];active.append(common.pop());edges.append(g.minus(c,a))
    assert len(set(active))==len(active) and len(active)>=3
    assert all(g.cross(l,r)>0 for l,r in zip(edges,edges[1:]))
    _,_,c3,b3,*_=g.setup()
    assert g.cross(g.rot(g.minus(v,c3)),edges[0])>0 and g.cross(edges[-1],g.rot(g.minus(t,b3)))>0
    assert all(g.cross(g.minus(t,v),g.minus(z,v))<0 for z,_ in verts[1:-1]);assert v[0]<t[0]
    pts=[z for z,_ in verts];loss=sum(g.cross(g.minus(a,v),g.minus(c,v)) for a,c in zip(pts[1:-1],pts[2:]))/2
    gain=g.arb(g.trial(params[-1][0],tight=True)['gain_area'])-loss-prev_gain;assert gain>0
    return dict(n=n,relaxed=relaxed,normals=normals,offsets=offsets,ys=ys,high=ys[0],v=v,t=t,u=u,gain=gain,active=active,breaks=[str(z[1]) for z in pts[1:-1]])


def sign_case(model,old,geo,params,key):
    targeted=key.startswith('U'); old_pair=key in old['ordered_branches']
    idx=int(key[1:]) if targeted else None
    depth='0.000555' if old_pair else (params[idx][0] if targeted else params[-1][0]); p=model(depth)
    basis=np.eye(p.dimension+1,dtype=object)
    letters=['U','1'] if targeted else list(key)
    for i,lbl in enumerate(letters):
        point=basis[2*i:2*i+2]
        if lbl in ('0','1'):
            j=int(lbl);p.linear(-np.asarray(geo['n'][j],dtype=object)@point,-geo['relaxed'][j])
            # In a mixed case this point really hit the old physical two-facet ear.
            if not old_pair:
                norm=b.unit_outer(geo['u'],geo['t']);p.linear(-np.asarray(norm,dtype=object)@point,-g.dot(norm,geo['t']))
                p.linear(point[1],geo['u'][1])
        else:
            height=geo['high'] if lbl=='H' else geo['ys'][idx]
            p.linear(-point[1],-height)
            if lbl=='U':
                j=idx+2;p.linear(-np.asarray(geo['normals'][j],dtype=object)@point,-geo['offsets'][j])
    return p


def run(archive,ah_file,param_file,weights,generate,out_weights):
    oldpath=archive/'eh-two-line-cover-certificate.json';old=json.loads(oldpath.read_bytes());params=json.loads(param_file.read_bytes());ah=json.loads(ah_file.read_bytes())
    assert old['parameters']==['0.000555','0.0005522','0.0075'] and ah['source_hashes'][oldpath.name]==sha256(oldpath.read_bytes()).hexdigest()
    model,model_hash=b.eh.specialize(F('1.5497'));assert model_hash==old['specialized_model_sha256'];geo=front(old,params)
    base=old['ordered_branches'];cases=[*base,'H0','0H','1H','HH',*(f'U{i}' for i in range(len(params)))];bounds={}
    for key in cases:
        p=sign_case(model,old,geo,params,key);initial=len(p.rows);obj=b.eh.phi_objective(p,b.RIGHT,1)
        if key in base or not generate:
            s=base[key] if key in base else weights['branches'][key];b.replay_cuts(p,s['adaptive_cuts']);du=s['dual']
        else:
            for ss in base.values():b.replay_cuts(p,ss['adaptive_cuts'])
            du=p.optimize(obj);s={'adaptive_cuts':list(p.adaptive_cuts),'dual':deepcopy(du)};b.trim(s,[s['dual']],initial);weights['branches'][key]=s
        upper=p.validate(obj,du)-g.arb(3)/2;assert upper<b.arbq(b.F_LIMIT),(key,str(upper));bounds[key]=str(upper.upper())
    p=model(params[-1][0]);b.eh.add_joint_query(p);p.linear(b.eh.phi_objective(p,b.RIGHT,1),g.arb(3)/2+b.arbq(b.F_LIMIT));initial=len(p.rows)
    stencil=deepcopy(old['enclosing_pose_certificate']) if generate else weights['pose'];b.replay_cuts(p,stencil['adaptive_cuts'])
    if generate:
        validate=p.validate
        def new_dual(obj,du):
            if du is not stencil['right_weights']:du[:]=p.optimize(obj)
            return validate(obj,du)
        p.validate=new_dual
    pose=b.eh.certify(p,b.RIGHT,10,stencil,eh_relocation=True,joint_eh=True,revisit_ah=True)
    if generate:
        duals=[pose['right_weights']]+[du for x in pose['edges'].values() for du in x['weights'].values()];b.trim(pose,duals,initial);pose.pop('summary');weights['pose']=pose;out_weights.write_text(json.dumps(weights,separators=(',',':'))+'\n');return
    u,v,t=(geo[k] for k in ('u','v','t'));cost=g.norm(g.minus(v,t))*g.norm(g.minus(b.unit_outer(v,t),b.unit_outer(u,t)))
    m=placement_lower(ah);assert m>cost
    clear=b.prev.ah_clearance(ah,[(a,c) for a in (v[0],t[0]) for c in (v[1],t[1])])
    gains=list(map(b.eh.gain_endpoint,(geo['gain'].lower(),geo['gain'].upper())));l,h=map(F,ah['area_exact_interval']);ends=(l-gains[1],h-gains[0])
    return {'AH_sha256':sha256(ah_file.read_bytes()).hexdigest(),'R126_sha256':sha256(oldpath.read_bytes()).hexdigest(),'parameter_sha256':sha256(param_file.read_bytes()).hexdigest(),'model_sha256':model_hash,
            'active_facet_indices':geo['active'],'interior_vertex_heights':geo['breaks'],'F_right_upper_by_case':bounds,'F_limit_exact':str(b.F_LIMIT),'EH_placement':pose['summary'],'arc_signs':b.eh.voronoi_signs(b.eh.arc_data()),
            'AH_local_chord_transfer_cost_upper':str(cost.upper()),'AH_after_transfer_lower':str((m-cost).lower()),'whole_rectangle_AH_clearance':clear,'EH_gain_over_R126_exact_interval':list(map(str,gains)),
            'gain_over_R126_arb':str(geo['gain']),'area_exact_interval':list(map(str,ends)),'area_decimal_interval':b.prev.format_area(*ends)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('AH',type=Path);p.add_argument('parameters',type=Path);p.add_argument('weights',type=Path);p.add_argument('output',type=Path);p.add_argument('--generate',action='store_true');a=p.parse_args()
    if a.generate:run(a.archive,a.AH,a.parameters,{'branches':{}},True,a.weights)
    result=run(a.archive,a.AH,a.parameters,json.loads(a.weights.read_bytes()),False,None);result['weights_sha256']=sha256(a.weights.read_bytes()).hexdigest();a.output.write_text(json.dumps(result,indent=2)+'\n')
