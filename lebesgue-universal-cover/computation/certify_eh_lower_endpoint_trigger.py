"""Certify a single constant-width body showing sharpness of the discrete F-prime trigger at T."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from certify_eh_joint_parameters import specialize
import audit_optimized_eh_chord as g


def verify():
    specialize(F('1.5497'));n,m,_,b3,_,_,s,t,_=g.setup();f=g.unit(g.PI/6-g.SIGMA+2*g.PI/3)
    r=g.line_pair(m[2],n[1],g.HALF,g.arb(0))
    pts={'S':s,'S-n0':g.minus(s,n[0]),'R':r,'R-f':g.minus(r,f),'B3':b3,'B3-n1':g.minus(b3,n[1]),
         'T':t,'n2/2':g.scale(g.HALF,n[2]),'-n2/2':g.scale(-g.HALF,n[2])}
    pairs=[('S','S-n0'),('S','R'),('S','T'),('R','R-f'),('B3','B3-n1'),('B3','T'),('n2/2','-n2/2')]
    exact={frozenset(v) for v in pairs};gaps=[]
    for (i,a),(j,b) in combinations(pts.items(),2):
        if frozenset((i,j)) in exact:continue
        d=g.minus(a,b);gaps.append((1-g.dot(d,d)).lower())
    assert min(gaps)>0
    radius=g.arb('0.001');other=min((1-g.norm(g.minus(t,q))-radius).lower() for k,q in pts.items() if k not in ('S','B3'))
    assert other>0 and 2*radius<1
    # Along the B3 short arc x=B3+R_angle(T-B3), |x-S|^2 has derivative
    # 2 (B3-S).J R_angle(T-B3). The continuation towards I is anticlockwise.
    theta=g.arb(0).union(g.arb('0.0009'));co,si=g.unit(theta);rad=g.minus(t,b3)
    radial=g.plus(g.scale(co,rad),g.scale(si,g.rot(rad)));derivative=2*g.dot(g.minus(b3,s),g.rot(radial))
    assert derivative<0
    # The next boundary endpoint I is the radial projection onto m3.x=1/2.
    from audit_eh_continuous_placement import project
    i=project(b3,m[3]);i_radial=g.minus(i,b3)
    assert g.cross(rad,i_radial)>0 and g.dot(rad,i_radial)>0 and g.cross(radial,i_radial)>0
    return {'sigma_degrees':'1.5497','lens_neighborhood_radius':'0.001','verified_arc_angle_interval':['0','0.0009'],
            'anchors':{k:list(map(str,x)) for k,x in pts.items()},'exact_unit_pairs':pairs,'other_pair_count':len(gaps),
            'minimum_other_pair_squared_slack_lower':str(min(gaps)),'minimum_other_anchor_distance_slack_after_radius_lower':str(other),
            'arc_squared_distance_derivative_upper':str(derivative.upper()),'arc_before_old_I_cross_lower':str(g.cross(radial,i_radial).lower()),
            'P_support_guards':{'n0':['S','S-n0'],'n1':['B3','B3-n1'],'n2':['n2/2','-n2/2'],'m3_negative':'B3','m5_negative':'R','Fprime_negative':'R-f','Cprime_negative':'R'},
            'precision_bits':g.ctx.prec}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args();a.output.write_text(json.dumps(verify(),indent=2)+'\n')
