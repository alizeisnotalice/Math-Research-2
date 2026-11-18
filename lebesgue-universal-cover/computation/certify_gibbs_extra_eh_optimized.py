"""Use the new U endpoint and signed projections to certify an extra EH cut."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

trial_delta = F(sys.argv[1]) if len(sys.argv)>1 else F(154,10**7)
assert trial_delta > 0
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
delta=trial_delta
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
    "C3_arc_enters_ear_side_at_U":-dot(u1t,sub(u1,c3)),
    "C3_arc_enters_ear_side_at_J":-dot(u1t,sub(j,c3)),
    "B3_arc_increases_to_chord_at_J":dot(u1t,sub(j,b3)),
    "B3_arc_increases_to_chord_at_T":dot(u1t,sub(t,b3)),
    "new_U_x_left_of_old_U":u[0]-u1[0],
    "new_U_y_higher_than_old_U":u1[1]-u[1],
    "new_U_shift_y_no_more_than_delta":I(delta)-(u1[1]-u[1]),
    "new_chord_ratio_x_to_minus_y_below_13_200":F(13,200)*(-u1t[1])-u1t[0],
    "new_chord_x_increasing":u1t[0],
    "new_chord_y_decreasing":-u1t[1],
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

# The new cut is below U1 and on the right of the chord U1->T.
# With y=U1_y-z_y>=0, z_x-U1_x<=k*y.  The two overlap points have
# total height loss <=2*eta/pal_slope by the unit diameter constraint.
k=u1t[0]/(-u1t[1])
eta=half-dot(p,u1)
pal_slope=half-r3*k/2
assert pal_slope.lo>0 and eta.lo>0
cone_checks={
    "positive_Pal_loss_at_U":eta,
    "Pal_height_slope":pal_slope,
}
cone_lower={key:strict(value,key) for key,value in cone_checks.items()}
mean_height=(eta/pal_slope).hi

# Use the new U1 on the unit sphere about C3, instead of the old U:
#   r(2*a+r)<=2*(horizontal*k-a-r)*y.
# Its increasing concave root f(y) bounds individual r's.  Jensen
# bounds their mean by f(mean_height).
a=c3[1]-u1[1]
horizontal=half-u1[0]
product=horizontal*k
assert a.lo>0 and (product-a).lo>0
def cb_shift_bound(y):
    return (a*a+2*product*y+y*y).sqrt()-a-y
r_bound=cb_shift_bound(I(mean_height)).hi
individual_r_bound=cb_shift_bound(I(2*mean_height)).hi
assert 0 < r_bound <= individual_r_bound

# e=b-w_K-m has parallel component in [-Mr,0]. More precisely,
# the tangent component is at most sqrt(2Mr-M^2r^2) when Mr<1.
# The signed loss in the chosen projection is bounded by
#   beta*sqrt(2Mr-M^2r^2) - lambda*r,
#   lambda=v_y + M*(v.m).
target=sub(c3,m[5])
normal,tangent=n[2],rot(n[2])
midpoint_b=mul(-half,normal)
relative=sub(target,midpoint_b)
m5_slope=-m[5][1]
rr=I(r_bound)*m5_slope
displacement=(2*rr-rr*rr).sqrt()
# The endpoint loss is linear in the measuring direction on the branch
# v.Jm>0, v.m<0. Align v with the corresponding corrected vector.
corrected=add(add(relative,(I(0),I(r_bound))),sub(mul(rr,m[5]),mul(displacement,rot(m[5]))))
corrected_length=norm2(corrected).sqrt()
test_direction=mul(1/corrected_length,corrected)
direction_tangent=dot(test_direction,rot(m[5]))
assert direction_tangent.lo>0
beta=direction_tangent
parallel_coefficient=dot(test_direction,m[5])
lambda_signed=test_direction[1]+m5_slope*parallel_coefficient
assert parallel_coefficient.hi<0 and m5_slope.lo>0
assert (test_direction[1]-pal_slope).lo>0 and lambda_signed.lo>0
assert (1-m5_slope*I(individual_r_bound)).lo>0
# The loss is concave and is increasing up to the certified individual r bound.
loss_derivative_lower=(
    beta*m5_slope*(1-m5_slope*I(individual_r_bound))/
    (2*I(individual_r_bound)*m5_slope).sqrt() - lambda_signed
)
assert loss_derivative_lower.lo>0
projected_error=(beta*displacement-lambda_signed*I(r_bound)).hi
overlap_margin=corrected_length-1
assert overlap_margin.lo>0,overlap_margin.out()

# Green integral over the additional cap C3 arc U1->U, old U->T chord
# and new T->U1 chord.  Its positive orientation is forced by convexity.
gain=(delta+cross(sub(c3,t),sub(u,u1)))/2
assert gain.lo>0, gain.out()
print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "coordinate_source_sha256":hashlib.sha256(raw).hexdigest(),
    "sigma_degrees_exact":str(sigma),
    "clockwise_U_shift_radians_exact":str(delta),
    "new_U_interval":[x.out() for x in u1],
    "geometry_strict_lower_bounds":geom_lower,
    "cone_strict_lower_bounds":cone_lower,
    "mean_overlap_height_below_new_U_upper":str(mean_height),
    "individual_CB_contact_displacement_upper":str(individual_r_bound),
    "loss_derivative_lower":str(loss_derivative_lower.lo),
    "signed_coefficient_lambda_lower":str(lambda_signed.lo),
    "projection_direction_y_minus_pal_slope_lower":str((test_direction[1]-pal_slope).lo),
    "mean_CB_contact_displacement_upper":str(r_bound),
    "mean_opposite_contact_projection_error_upper":str(projected_error),
    "projected_distance_overlap_contradiction_lower":str(overlap_margin.lo),
    "additional_cut_area_interval":gain.out(),
    "passed":True,
},indent=2))
