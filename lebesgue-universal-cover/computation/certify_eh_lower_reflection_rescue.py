"""Simultaneous failure of identity/Q1 at the lower tip and a uniform 300-degree rescue for this family."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
from certify_eh_all_completions import Lens
from certify_eh_general_relocation import rotate_angle
from certify_eh_joint_parameters import specialize
from certify_eh_polyline_obstruction import clip_geometry
import audit_optimized_eh_chord as g


def verify():
    specialize(F('1.5497'));n,m,_,b3,_,_,s,t,_=g.setup();f=g.unit(g.PI/6-g.SIGMA+2*g.PI/3);r=g.line_pair(m[2],n[1],g.HALF,g.arb(0))
    def reflect(z):return rotate_angle((z[0],-z[1]),5*g.PI/3)
    def inv(z):return rotate_angle(z,-5*g.PI/3)
    pts={'S':s,'E':g.minus(s,n[0]),'R':r,'R-f':g.minus(r,f),'B':b3,'B-n1':g.minus(b3,n[1]),'T':t}
    pts.update({f'Q({k})':reflect(pts[k]) for k in ('S','E','R-f','T')})
    unit=[('S','E'),('S','R'),('S','T'),('R','R-f'),('B','B-n1'),('B','T'),('R','Q(R-f)'),('R','Q(S)'),('Q(S)','Q(T)'),('Q(S)','Q(E)'),('B-n1','Q(T)')]
    assigned={frozenset(q) for q in unit};gaps=[(1-g.dot(g.minus(a,b),g.minus(a,b))).lower() for (i,a),(j,b) in combinations(pts.items(),2) if frozenset((i,j)) not in assigned]
    assert min(gaps)>0
    radius=g.arb('.001');lower=min((1-g.norm(g.minus(t,z))-radius).lower() for k,z in pts.items() if k not in ('S','B'))
    separation=g.norm(g.minus(t,reflect(t)))+2*radius;assert lower>0 and separation<1
    v=clip_geometry(['.00055845','.0005522','.0075'])[3];L=Lens(pts)
    direction=g.minus(b3,t);det=g.cross(n[1],direction);mix=[g.cross(m[0],direction)/det,g.cross(n[1],m[0])/det];assert min(mix)>0
    positive_cap=g.HALF-L.upper_support(inv(m[3]));cprime=g.HALF-L.upper_support(inv(g.neg(f)));assert min(positive_cap,cprime)>0
    phi=g.arb('.0009');a=g.arb(0).union(phi);rr=g.minus(t,b3);vr=rotate_angle(rr,a);p=g.plus(b3,vr)
    det=g.cross(g.minus(p,v),g.rot(vr));assert det>0
    d0=g.minus(t,v);dend=g.minus(g.plus(b3,rotate_angle(rr,phi)),v);forward=g.dot(d0,g.minus(p,v));assert forward>0 and g.cross(d0,dend)>0
    margins=[]
    for d in [d0,dend]:
        normal=g.scale(-1/g.norm(d),g.rot(d));margins.append((g.dot(normal,v)-L.upper_support(inv(normal))).lower())
    assert min(margins)>0
    return {'sigma_degrees':'1.5497','EH_upper_depth':'.00055845','lower_arc_angle_range':['0','.0009'],'local_radius':'.001','anchors':{k:list(map(str,z)) for k,z in pts.items()},
            'unit_pairs':unit,'remaining_pairs':len(gaps),'remaining_squared_distance_slack_lower':str(min(gaps)),'other_anchor_slack_after_radius':str(lower),
            'two_neighborhoods_max_distance':str(separation.upper()),'m0_positive_cone_at_B3':list(map(str,mix)), 'rescue_rotation_degrees':300,
            'rescue_positive_E_slack_lower':str(positive_cap.lower()),'rescue_negative_C_slack_lower':str(cprime.lower()),'rotating_chord_angle_derivative_numerator_lower':str(det.lower()),
            'endpoint_chord_slacks_lower':list(map(str,margins))}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);s=p.parse_args();s.output.write_text(json.dumps(verify(),indent=2)+'\n')
