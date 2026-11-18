"""Move both EH support faces by the available ordered-contact dual budgets."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from flint import __version__ as flint_version
import certify_eh_two_line_cover as eh
from certify_ah_three_line_cut import format_area
from certify_eh_general_relocation import exact_fraction

g = eh.g
RELAX = ('0.00000001036', '0.0000000000086')


def point_on_faces(normals, offsets):
    a, b = normals
    h, k = offsets
    det = g.cross(a, b)
    return ((h*b[1] - k*a[1])/det, (a[0]*k-b[0]*h)/det)


def ah_clearance(ah, points):
    from audit_eh_continuous_placement import circle_pair, line_pair
    n,m,*_ = g.setup()
    a = line_pair(n[4],m[3]); c=(-a[0],a[1]); root=g.arb(3).sqrt()
    l0=(root/2-g.HALF,(root/2-g.HALF-1)/root)
    z=circle_pair(l0,(-g.HALF/2,-root/4),-1)
    v=g.minus(z,c); v2=g.dot(v,v)
    ta=g.plus(c,g.scale(1/v2,g.minus(v,g.scale((v2-1).sqrt(),g.rot(v)))))
    D=g.minus(ta,z); axis=g.scale(-1/g.norm(D),D); normal=g.rot(axis)
    cuts=[('0.355','0.00060'),('0.410','0.00030'),('0.404','0.00046')]+list(ah['lines'])
    slacks=[]
    for p,angle in cuts:
        c,s=g.unit(exact_fraction(Fraction(angle))); line=g.plus(g.scale(c,axis),g.scale(s,normal))
        pivot=g.plus(z,g.scale(exact_fraction(Fraction(p)),D))
        slacks.extend(g.cross(line,g.minus(x,pivot)) for x in points)
    lower=min(v.lower() for v in slacks)
    assert lower>0
    return str(lower)


def verify(archive, ah_path):
    original = archive/'eh-two-line-cover-certificate.json'
    old = json.loads(original.read_bytes())
    ah = json.loads(ah_path.read_bytes())
    assert ah['R167_area_sha256']==sha256((archive/'ah-three-line-cut-certificate.json').read_bytes()).hexdigest()
    assert old['parameters']==['0.000555','0.0005522','0.0075']
    model, model_hash = eh.specialize(Fraction('1.5497'))
    assert model_hash == old['specialized_model_sha256']
    right = Fraction(old['right_degrees'])
    normals, offsets, _, u, w = eh.clip_geometry(old['parameters'])
    t = g.setup()[7]
    relax = [exact_fraction(Fraction(s)) for s in RELAX]
    reduced = [b-r for b,r in zip(offsets,relax)]
    x = point_on_faces(normals,reduced)
    # The old and new hinges are strictly on the same side of the enclosing chord.
    # Positive barycentric coordinates certify x belongs to the original triangle.
    uw, wt = g.minus(w,u),g.minus(t,w)
    turn = g.cross(uw,wt)
    distances = [g.cross(g.minus(b,a),g.minus(x,a)) for a,b in ((u,w),(w,t),(t,u))]
    assert turn>0 and all(value>0 for value in distances)
    gain = (distances[0]+distances[1])/2
    assert gain>0
    directions = [g.minus(x,u),g.minus(t,x)]
    det = g.cross(*normals)
    mixture = []
    for v in directions:
        outer=g.neg(g.rot(v))
        coeff=(g.cross(outer,normals[1])/det,g.cross(normals[0],outer)/det)
        assert min(coeff)>0
        mixture.append([str(c) for c in coeff])
    _, _, c3, b3, *_ = g.setup()
    assert g.cross(g.rot(g.minus(u,c3)),directions[0])>0
    assert g.cross(directions[1],g.rot(g.minus(t,b3)))>0
    assert g.cross(*directions)>0 and all(v[1]<0 for v in directions)
    bounds=[]
    branch={}
    for key, seed in old['ordered_branches'].items():
        program=model(old['parameters'][0]); basis=np.eye(program.dimension+1,dtype=object)
        for i,j in enumerate(map(int,key)):
            support=np.array(normals[j],dtype=object)@basis[2*i:2*i+2]
            program.linear(-support,-reduced[j])
        for index,direction in seed['adaptive_cuts']:
            program.cut(index,direction,record=True)
        objective=eh.phi_objective(program,right,1)
        upper=program.validate(objective,seed['dual'])
        f=upper-g.arb(3)/2
        assert f<0,(key,str(f))
        bounds.append(upper); branch[key]=str(f.upper())
    program=model(old['parameters'][0]);eh.add_joint_query(program)
    program.linear(eh.phi_objective(program,right,1),max(bounds))
    seed=old['enclosing_pose_certificate']
    for index,direction in seed['adaptive_cuts']:
        program.cut(index,direction,record=True)
    placement=eh.certify(program,right,10,seed,eh_relocation=True,joint_eh=True,revisit_ah=True)
    low,high=map(Fraction,ah['area_exact_interval'])
    gains=[eh.gain_endpoint(v) for v in (gain.lower(),gain.upper())]
    endpoints=[low-gains[1],high-gains[0]]
    return {'base_EH_sha256':sha256(original.read_bytes()).hexdigest(),'base_AH_sha256':sha256(ah_path.read_bytes()).hexdigest(),
            'model_sha256':model_hash,'flint_version':flint_version,'precision_bits':g.ctx.prec,'old_parameters':old['parameters'],
            'unit_normal_relaxations_exact':list(RELAX),'new_hinge':list(map(str,x)), 'positive_coefficients':mixture,
            'minimum_AH_clearance_on_whole_old_EH_triangle':ah_clearance(ah,(u,w,t)),
            'triangle_edge_determinants':list(map(str,distances)), 'upper_F_by_ordered_branch':branch,
            'placement':placement['summary'],'area_gain_exact_interval':list(map(str,gains)),
            'area_gain_interval_arb':str(gain),'area_exact_interval':list(map(str,endpoints)),
            'area_decimal_interval':format_area(*endpoints)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path); parser.add_argument('baseline_AH',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();result=verify(args.archive,args.baseline_AH)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
