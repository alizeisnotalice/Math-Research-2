"""Certify the common EH relocation for a separate cap on the remaining Tp-P arc."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import audit_optimized_eh_chord as g
from audit_eh_continuous_placement import polygon
from certify_eh_general_relocation import rotate_angle
from certify_eh_joint_parameters import specialize


def verify():
    specialize(Fraction('1.5497'))
    n,m,*_=g.setup()
    f3=g.line_pair(n[3],m[2])
    center=g.circle_line(f3,n[5],1)
    d1=g.line_pair(n[4],n[5])
    e=g.minus(d1,center)
    e=g.scale(1/g.norm(e),e)
    _,bounds=polygon('.000555',1)
    tp_rad=next(v for k,v,h in bounds if k=='tangent_Tp_0')
    # At P the radius is n2, and the arc approaches P counterclockwise.
    arc=g.arb.atan2(g.cross(tp_rad,n[2]),g.dot(tp_rad,n[2]))
    base=g.arb('.022')
    half=g.arb('.0085')
    assert 0<base-half<base+half<arc
    sag=1-half.cos()
    # The perpendicular to the radius n2 points from G to D1 clockwise or
    # counterclockwise. Use its sign to evaluate the radial drift exactly.
    assert g.dot(e,g.rot(n[2]))>0
    # e = J n2 by the exact n5 face identity, so -e.R_(-t)n2 = sin(t).
    drift=(base-half).sin()
    tau=(2*sag+sag**2)/(2*drift)
    psi=g.arb(0).union(2*half)+base-half
    radial=rotate_angle(n[2],-psi)
    arc_point=g.plus(center,radial)
    mp=[g.unit(g.PI/6-g.SIGMA+j*g.PI/3) for j in range(6)]
    zguard=g.dot(mp[1],arc_point)-g.HALF-sag
    aguard=g.dot(mp[5],center)-g.HALF-tau
    reflection=lambda p:rotate_angle((p[0],-p[1]),2*g.PI/3)
    distance=g.norm(g.minus(center,reflection(center)))-2*tau-1
    gain=(2*half-(2*half).sin())/2
    assert min(zguard,aguard,distance,gain)>0
    return {'sigma_degrees':'1.5497','precision_bits':g.ctx.prec,'clockwise_center_from_P_radians':'.022',
            'angular_halfwidth_radians':'.0085','true_full_Tp_P_sweep':str(arc),
            'radial_sagitta_upper':str(sag.upper()),'G_contact_distance_upper':str(tau.upper()),
            'prime_1_hit_margin_lower':str(zguard.lower()),'prime_5_hit_margin_lower':str(aguard.lower()),
            'double_new_hit_diameter_excess_lower':str(distance.lower()),'true_circular_cap_area':str(gain),
            'remaining_case_for_full_AH_cover':'original hits this new arc cap, QA image hits prior AH'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(verify(),indent=2)+'\n')
