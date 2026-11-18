"""Prove that the entire Sprague arc T-I is mandatory inside the P host."""
import argparse
import json
from pathlib import Path
import audit_optimized_eh_chord as g
from certify_eh_general_relocation import rotate_angle
from certify_eh_mandatory_lower_lens import verify as verify_lens


def verify():
    prior=verify_lens()
    radius=g.arb('.05')
    slack=g.arb(prior['other_anchor_margin_after_radius_lower'])+g.arb('.001')-radius
    assert slack>0
    n,m,_,b,_,_,s,t,_=g.setup()
    v=g.minus(t,b)
    w=g.neg(m[0])  # I = B3 - m0
    assert g.cross(v,w)>0 and g.dot(v,w)>0
    angle=g.arb.atan2(g.cross(v,w),g.dot(v,w))
    max_distance=g.norm(g.minus(w,v))
    assert 0<max_distance<angle<radius
    sweep=g.arb(0).union(angle)
    derivative=2*g.dot(g.minus(b,s),g.rot(rotate_angle(v,sweep)))
    assert derivative<0
    return {
        'verified_radius':'.05',
        'other_anchor_distance_slack_on_ball_lower':str(slack.lower()),
        'T_I_sweep_radians':str(angle),
        'T_I_maximum_distance_from_T':str(max_distance),
        'squared_distance_to_S_derivative_upper':str(derivative.upper())
    }


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(verify(),indent=2)+'\n')
