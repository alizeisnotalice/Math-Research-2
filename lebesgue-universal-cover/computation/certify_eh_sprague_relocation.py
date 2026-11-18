#!/usr/bin/env python3
"""Certify EH relocation with all five inherited Sprague arc supports retained exactly."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

import audit_optimized_eh_chord as geometry
import certify_eh_relocation_adaptive as adaptive
from audit_eh_continuous_placement import polygon
from certify_eh_general_relocation import numeric
from certify_eh_joint_parameters import gain_endpoint, specialize


def arc_data():
    n,m,c3,b3,_,_,_,_,j=geometry.setup()
    meet=geometry.line_pair
    b1,d1=meet(n[0],n[1]),meet(n[4],n[5])
    b2,e3=meet(n[0],m[0]),meet(n[4],m[3])
    f3,f2=meet(n[3],m[2]),meet(n[2],m[2])
    g=geometry.circle_line(f3,n[5],1)
    x=geometry.circle_pair(g,e3,-1)
    return (('CS',f3,f2,m[5],g),('ES-B',b3,b2,m[3],j),
            ('ES-C',c3,b1,n[3],j),('AS-G',g,d1,n[2],x),('AS-E',e3,d1,n[1],x))


def voronoi_signs(data):
    dot,minus=geometry.dot,geometry.minus
    result={}
    for name,center,other,tangent_radial,strict_point in data:
        distant_radial=minus(strict_point,center)
        strict_loss=-dot(distant_radial,minus(other,center))
        minor_arc=dot(tangent_radial,distant_radial)
        assert strict_loss>0 and minor_arc>0
        # The other endpoint radial is the projection direction at the
        # opposite support line, orthogonal to the contact segment exactly.
        result[name]={'strict_endpoint_Voronoi_margin':str(strict_loss.lower()),
                      'endpoint_radial_dot_lower':str(minor_arc.lower())}
    return result


def add_contacts(program,data):
    previous=program.dimension
    program.dimension+=4
    eye=np.eye(program.dimension+1,dtype=object)

    def expand(point):
        return np.column_stack((point[:,:previous],np.zeros((2,4),dtype=object),point[:,-1]))

    program.witnesses=[expand(p) for p in program.witnesses]
    program.support=[expand(p) for p in program.support]
    program.disks=[expand(p) for p in program.disks]
    program.rows=[(np.concatenate((a,np.zeros(4,dtype=object))),b) for a,b in program.rows]
    program.numeric_rows=[numeric(a) for a,_ in program.rows]
    n,m,*_=geometry.setup()
    axis=n[5]
    reflection=2*np.outer(axis,axis)-np.eye(2,dtype=object)
    old_disks=len(program.disks)
    contact_by_name={name:(a,b) for name,a,b,_,_ in data}
    g,d1=contact_by_name['AS-G']
    direction=np.asarray(geometry.minus(g,d1),dtype=object)
    for q in (program.witnesses[5],reflection@program.witnesses[5]):
        shifted=q-np.outer(d1,eye[-1])
        program.linear(-direction@shifted,0)
        program.linear(direction@shifted,geometry.dot(direction,direction))

    for i,name in enumerate(('CS','ES-B')):
        a,b=contact_by_name[name]
        displacement=geometry.minus(b,a)
        for side in range(2):
            s=eye[previous+2*i+side]
            program.linear(-s,0)
            program.linear(s,1)
            point=np.outer(a,eye[-1])+np.outer(displacement,s)
            if side:
                point=reflection@point
            for normal in n+[m[3],m[5],reflection@m[3],reflection@m[5]]:
                program.linear(np.asarray(normal,dtype=object)@point,geometry.HALF)
            program.disks.extend(point-p for p in program.witnesses)
            program.witnesses.append(point)

    program.disk_num=np.stack([numeric(d) for d in program.disks])
    for index in range(old_disks,len(program.disks)):
        for direction in ((1,0),(-1,0),(0,1),(0,-1)):
            program.cut(index,direction)


PREFIXES=('chord_Q_','chord_Tp_','chord_H_','chord_Te_','chord_G_')


def hybrid_polygon(delta,subdivisions):
    inner,outer=polygon(delta,subdivisions)
    return [edge for edge in inner if not edge[0].startswith(PREFIXES)],outer


def run(args):
    sigma=Fraction(args.sigma)
    model_type,model_hash=specialize(sigma)
    checks=geometry.trial(args.delta,tight=True)
    data=arc_data()
    signs=voronoi_signs(data)
    program=model_type(args.delta)
    if args.contacts:
        add_contacts(program,data)
    adaptive.polygon=hybrid_polygon
    defaults=adaptive.multipliers

    def multipliers(name):
        if name=='K-C3':
            return tuple(map(Fraction,('0','1.155','1.35','1.5','2','2.5','3','5','8','20')))
        return defaults(name)

    adaptive.multipliers=multipliers
    replay=None
    if args.replay:
        replay=json.loads(args.replay.read_text())
        receipt=replay['summary']
        assert Fraction(receipt['sigma_degrees'])==sigma and receipt['delta']==args.delta
        assert Fraction(receipt['angle_root_range_degrees'][1])==Fraction(args.right)
        assert bool(receipt.get('explicit_Sprague_contacts',False))==args.contacts
        assert receipt['specialized_model_sha256']==model_hash
        for index,direction in replay['adaptive_cuts']:
            program.cut(index,direction,record=True)
    result=adaptive.certify(program,Fraction(args.right),args.max_depth,replay)
    summary=result['summary']
    summary['sigma_degrees']=str(sigma)
    summary['specialized_model_sha256']=model_hash
    summary['Sprague_arc_endpoint_signs']=signs
    summary['explicit_Sprague_contacts']=args.contacts
    gain=geometry.arb(checks['gain_area'])
    summary['cut_gain_area_interval_exact']=[str(gain_endpoint(gain.lower())),str(gain_endpoint(gain.upper()))]
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sigma')
    parser.add_argument('delta')
    parser.add_argument('right')
    parser.add_argument('output',type=Path)
    parser.add_argument('--contacts',action='store_true')
    parser.add_argument('--max-depth',type=int,default=10)
    parser.add_argument('--replay',type=Path)
    run(parser.parse_args())
