"""Locate the entire new lower EH ear near T and measure its exact true-arc gain."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import certify_eh_upper_chain as chain
from certify_eh_general_relocation import rotate_angle

g=chain.g


def verify(archive,param_file,ah_file):
    chain.b.eh.specialize(F('1.5497'));old=json.loads((archive/'eh-two-line-cover-certificate.json').read_bytes());params=json.loads(param_file.read_bytes());front=chain.front(old,params)
    _,_,_,b3,_,_,_,t,_=g.setup();first_n,first_offsets,_,_,_=chain.b.eh.clip_geometry(old['parameters']);hr=[h-chain.b.arbq(F(s)) for h,s in zip(first_offsets,chain.b.prev.RELAX)];x=chain.b.prev.point_on_faces(first_n,hr)
    alpha=chain.b.arbq(F('0.0000001'));lam=chain.b.arbq(F('0.00004'));y=g.plus(t,g.scale(lam,g.minus(x,t)));d=g.plus(b3,rotate_angle(g.minus(t,b3),alpha));new_edge=g.minus(d,y)
    assert g.cross(g.minus(t,x),new_edge)>0 and g.cross(new_edge,g.rot(g.minus(d,b3)))>0
    assert g.cross(new_edge,g.minus(t,y))<0 and g.cross(g.minus(t,x),g.minus(d,t))>0
    # At all old EH vertices above Y, the new supporting line is redundant
    # relative to the old lower supporting line since the two meet at Y.
    assert x[1]>y[1]>t[1]>d[1] and g.cross(new_edge,g.minus(x,y))>0
    prior_slacks=[h-g.dot(n,d) for n,h in zip(front['normals'],front['offsets'])];assert min(prior_slacks)>0
    cap_radius=max(g.norm(g.minus(y,t)).upper(),alpha.upper());assert cap_radius<g.arb('0.0000005')
    triangle=g.cross(g.minus(t,y),g.minus(d,y))/2;segment=(alpha-alpha.sin())/2;gain=triangle+segment
    assert triangle>0 and segment>0
    mix=chain.b.arbq(F('0.00000001'))
    mixed_tip=g.scale(1-mix,t)
    mixed_violation=-g.cross(new_edge,g.minus(mixed_tip,y))/g.norm(new_edge)
    assert mixed_violation>0
    ah=json.loads(ah_file.read_bytes());rect=[(a,b) for a in (front['v'][0],d[0]) for b in (d[1],front['v'][1])];clear=chain.b.prev.ah_clearance(ah,rect)
    return {'explicit_mixing_weight':str(F('0.00000001')),'mixed_tip_distance_outside_new_bevel_lower':str(mixed_violation.lower()),'lower_depth':str(F('0.0000001')),'hinge_fraction':str(F('0.00004')),'radius_bound':str(cap_radius),'triangle_area':str(triangle),'true_arc_segment_area':str(segment),'true_area_gain':str(gain),'area_gain_exact_interval':[str(chain.b.eh.gain_endpoint(q)) for q in (gain.lower(),gain.upper())],
            'minimum_old_facet_slack_at_new_lower_endpoint':str(min(s.lower() for s in prior_slacks)),'AH_clearance_on_full_rectangle':clear,'AH_sha256':sha256(ah_file.read_bytes()).hexdigest(),'parameters_sha256':sha256(param_file.read_bytes()).hexdigest()}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('parameters',type=Path);p.add_argument('AH',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.write_text(json.dumps(verify(a.archive,a.parameters,a.AH),indent=2)+'\n')
