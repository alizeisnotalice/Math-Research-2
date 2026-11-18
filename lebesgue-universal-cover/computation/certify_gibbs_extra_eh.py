"""Exact signs, overlap exclusion and area gain for one larger EH cut."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

source = Path(__file__).with_name("certify_gibbs_area_refined.py")
raw = source.read_bytes()
parts = raw.decode().split("width = F(1,cells)")
assert len(parts)==2
sys.argv = [str(source), "1.5494009336023418793", "2", "1"]
d = {"__file__":str(source)}
exec(compile(parts[0],str(source),"exec"),d)
I,add,sub,mul,dot,cross,rot,cc = (d[k] for k in ("I","add","sub","mul","dot","cross","perp","cc"))
half,r3 = d["half"],d["r3"]
n,m = d["n"],d["m"]
sigma = d["sigma"]
trig,pi=d["trig"],d["pi_machin"]
u,t,s,c3,b3,h,ip=(d[k] for k in ("u","t_eh","s","c3","b3","h","ip"))
j=cc(b3,c3,-1)
delta=F(3,10**8)
co,si=(I(*trig.trig((delta,delta),kind)) for kind in ("cos","sin"))
rad=sub(u,c3)
u1=add(c3,sub(mul(co,rad),mul(si,rot(rad))))
p=(-r3/2,I(-half))
arg_e=trig.scale(pi,F(7,6)-sigma/180)
e_neg=tuple(I(*trig.trig(arg_e,kind)) for kind in ("cos","sin"))

def norm2(x):
    return dot(x,x)

def strict(x,label):
    assert x.lo>0,(label,x.out())
    return str(x.lo)

def minor(a,b,c):
    va,vb=sub(a,c),sub(b,c)
    assert cross(va,vb).lo>0 and dot(va,vb).lo>0

minor(h,u1,c3)
minor(u1,u,c3)
minor(u1,j,c3)
minor(j,t,b3)
u1t=sub(t,u1)
geom={
    "new_U_left_turn":cross(rot(sub(u1,c3)),u1t),
    "new_T_left_turn":cross(u1t,rot(sub(t,b3))),
    "new_U_forward_join":dot(rot(sub(u1,c3)),u1t),
    "new_T_forward_join":dot(u1t,rot(sub(t,b3))),
    "S_disk_departure_from_T":dot(sub(t,s),sub(u1,t)),
    "C3_arc_S_disk_increasing_at_new_U":2*dot(sub(c3,s),rot(sub(u1,c3))),
    "C3_arc_S_disk_increasing_at_J":2*dot(sub(c3,s),rot(sub(j,c3))),
    "B3_arc_S_disk_decreasing_at_J":-2*dot(sub(b3,s),rot(sub(j,b3))),
    "B3_arc_S_disk_decreasing_at_T":-2*dot(sub(b3,s),rot(sub(t,b3))),
    "new_U_outside_S_disk":norm2(sub(u1,s))-1,
    "Eprime_U_slack":dot(e_neg,u1)-half,
    "Eprime_T_slack":dot(e_neg,t)-half,
    "Eprime_C3_derivative_U":dot(e_neg,rot(sub(u1,c3))),
    "Eprime_C3_derivative_J":dot(e_neg,rot(sub(j,c3))),
    "negative_Eprime_B3_derivative_J":-dot(e_neg,rot(sub(j,b3))),
    "negative_Eprime_B3_derivative_T":-dot(e_neg,rot(sub(t,b3))),
    "new_U_y_below_S":s[1]-u1[1],
    "T_y_below_S":s[1]-t[1],
    "C3_arc_y_decreasing_at_U":-(u1[0]-c3[0]),
    "C3_arc_y_decreasing_at_J":-(j[0]-c3[0]),
    "B3_arc_y_decreasing_at_J":-(j[0]-b3[0]),
    "B3_arc_y_decreasing_at_T":-(t[0]-b3[0]),
    "S_outside_new_EH_by_x_bound":s[0]-t[0],
    "C3_arc_x_increasing_at_U":-(u1[1]-c3[1]),
    "C3_arc_x_increasing_at_J":-(j[1]-c3[1]),
    "B3_arc_x_increasing_at_J":-(j[1]-b3[1]),
    "B3_arc_x_increasing_at_T":-(t[1]-b3[1]),
}
geom_lower={k:strict(v,k) for k,v in geom.items()}

# All of both circle arcs and the chord lie in a cone of opening
# cos^{-1}(2/5) about the Pál normal from the new U.
cone=F(2,5)
eta=half-dot(p,u1)
arc_box=(I(j[0].lo,t[0].hi),I(t[1].lo,j[1].hi))
box_dist=norm2(sub(arc_box,u1)).hi
palt=dot(p,u1t)
cone_checks={
    "positive_Pal_loss_at_U":eta,
    "Pal_loss_no_more_than_half_delta":I(delta/2)-eta,
    "C3_arc_Pal_derivative_U":dot(p,rot(sub(u1,c3)))-cone,
    "C3_arc_Pal_derivative_J":dot(p,rot(sub(j,c3)))-cone,
    "B3_arc_Pal_derivative_J_negative":-dot(p,rot(sub(j,b3))),
    "B3_arc_Pal_derivative_T_negative":-dot(p,rot(sub(t,b3))),
    "Pal_chord_advance":palt,
    "Pal_chord_cone":palt*palt-cone*cone*norm2(u1t),
    "Pal_B3_arc_cone":palt*palt-cone*cone*I(box_dist),
}
cone_lower={k:strict(v,k) for k,v in cone_checks.items()}

# A simultaneous Q2 hit forces both hit locations within E=3.5*delta of U.
# A CB contact b then differs from C3 by at most r_bound in the y coordinate.
# The forced opposite m5 contact lies within rho of C3-m5.
epsilon=F(7,2)*delta
a=c3[1]-u[1]
assert a.lo>0 and (-m[5][1]).lo>0
r_bound=(I(2*epsilon+epsilon*epsilon)/(2*a)).hi
rho=(I(r_bound)+(2*I(r_bound)*(-m[5][1])).sqrt()).hi
target=sub(c3,m[5])
normal,tangent=n[2],rot(n[2])
parallel=dot(normal,target)
transverse=dot(tangent,target)
lower_distance2=(parallel+half)*(parallel+half)+transverse*transverse
overlap_margin=lower_distance2-(1+rho)**2
assert overlap_margin.lo>0,overlap_margin.out()

# Green integral over the additional cap C3 arc U1->U, old U->T chord
# and new T->U1 chord.  Its positive orientation is forced by convexity.
gain=(delta+cross(sub(c3,t),sub(u,u1)))/2
assert gain.lo>F(6,10**12), gain.out()
print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "coordinate_source_sha256":hashlib.sha256(raw).hexdigest(),
    "sigma_degrees_exact":str(sigma),
    "clockwise_U_shift_radians_exact":str(delta),
    "new_U_interval":[x.out() for x in u1],
    "geometry_strict_lower_bounds":geom_lower,
    "cone_strict_lower_bounds":cone_lower,
    "overlap_points_distance_from_original_U_upper":str(epsilon),
    "CB_contact_displacement_upper":str(r_bound),
    "opposite_contact_displacement_upper":str(rho),
    "squared_distance_overlap_contradiction_lower":str(overlap_margin.lo),
    "additional_cut_area_interval":gain.out(),
    "passed":True,
},indent=2))
