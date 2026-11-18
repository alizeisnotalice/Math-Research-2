"""Quantitative rescue of three cap-sign classes for arbitrary constant-width bodies near T."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_general_relocation import rotate_angle
from certify_eh_joint_parameters import specialize
from certify_eh_polyline_obstruction import clip_geometry


def outer(a,b):
    d=g.minus(b,a);return g.scale(-1/g.norm(d),g.rot(d))


def verify(ah_file):
    specialize(F('1.5497'));n,m,c3,b3,p,_,s,t,j=g.setup();v=clip_geometry(['.00055845','.0005522','.0075'])[3];u=clip_geometry(['.000555','.0005522','.0075'])[3]
    r=g.line_pair(m[2],n[1],g.HALF,g.arb(0));f3=g.line_pair(n[3],m[2]);b2=g.line_pair(n[0],m[0]);e=g.minus(b2,b3);e=g.scale(1/g.norm(e),e)
    assert r[1]>s[1]>c3[1] and all(g.dot(g.minus(b,r),g.minus(f3,r))<0 for b in (c3,s))
    radius=g.arb('0.0000005');D=2*radius+radius**2;cb=D/(2*(s[1]-t[1]));tau=D/(2*g.dot(g.minus(b3,t),e));eta=tau+(2*tau*abs(g.dot(n[1],e))).sqrt()
    assert 0<cb<g.arb('.000016') and 0<tau<g.arb('.000037') and 0<eta<g.arb('.006')
    eprime=g.unit(7*g.PI/6-g.SIGMA);cap=g.dot(eprime,t)-g.HALF-radius;assert cap>0
    theta=g.arb('.0009');a=g.arb(0).union(theta);zz=g.minus(t,b3);rad=rotate_angle(zz,a);pt=g.plus(b3,rad)
    end=g.plus(b3,rotate_angle(zz,theta));da=g.minus(pt,v);rotation=g.cross(da,g.rot(rad));assert rotation>0 and g.dot(g.minus(t,v),da)>0
    assert rad[0]<0 and rad[1]<0 and v[1]>t[1]>end[1] and v[0]<t[0]<end[0]
    # This coarse rectangle contains the full arc cap for every lower endpoint.
    assert p[0]<0 and p[1]<0
    pal=g.dot(p,t)-radius+g.dot(p,(end[0],v[1]))-1;assert pal>0
    slacks={}
    for label,te in [('zero',t),('maximum',end)]:
        nn=outer(v,te);h=g.dot(nn,v)
        reflection=h-1-g.dot((nn[0],-nn[1]),s)-cb
        turn=h-1-g.dot(rotate_angle(nn,-2*g.PI/3),g.minus(b3,n[1]))-eta
        assert reflection>0 and turn>0
        slacks[label]={'horizontal_reflection_lower':str(reflection.lower()),'rotation120_lower':str(turn.lower())}
    # To carry the AH common root to a changing lower endpoint it is enough
    # to estimate each change of unit normal on the cap of a putative violator.
    small=g.arb('.00002');L=g.norm(g.minus(t,v));length=L+small
    cost=g.norm(g.minus(outer(u,t),outer(v,t)))*(length+small)+2*small*length/(L-small)
    ah=json.loads(ah_file.read_bytes());rows=ah['bins']+ah.get('root_conditioned_classes',[])
    ids=[row['bin'] if 'bin' in row else row['reference_bin'] for row in rows];assert sorted(ids)==list(range(68))
    margin=min(g.arb(x['placement_margin_lower']).lower() for x in rows);assert margin>cost
    return {'radius':'.0000005','upper_EH_depth':'.00055845','old_AH_lower_arc_interval':['0','.0009'],'full_AH_lower_arc_interval':['0','.00002'],
            'anchor_S_distance_upper':str(cb.upper()),'B3_segment_parameter_upper':str(tau.upper()),'B3_minus_n1_distance_upper':str(eta.upper()),
            'Eprime_hit_margin_lower':str(cap.lower()),'Q2_pal_width_excess_lower':str(pal.lower()),'two_endpoint_support_slacks':slacks,
            'chord_angle_derivative_numerator_lower':str(rotation.lower()),'AH_sha256':sha256(ah_file.read_bytes()).hexdigest(),'AH_reference_classes':len(rows),'AH_root_transfer_cost_upper':str(cost.upper()),'AH_after_transfer_slack_lower':str((margin-cost).lower())}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('AH',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.write_text(json.dumps(verify(a.AH),indent=2)+'\n')
