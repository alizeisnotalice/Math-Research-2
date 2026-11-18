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
delta=F(11,10**6)
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

# The new cut is below U1 and lies on the right of the oriented
# chord U1->T, so for d_y=U1_y-z_y>=0, we have z_x-U1_x <= k*d_y.
# The Pál projection thus controls the sum of the two d_y's under overlap.
k=F(13,200)
eta=half-dot(p,u1)
pal_slope=half-r3*k/2
assert pal_slope.lo>0
cone_checks={
    "positive_Pal_loss_at_U":eta,
    "Pal_loss_no_more_than_half_delta":I(delta/2)-eta,
    "Pal_height_slope":pal_slope,
}
cone_lower={k:strict(v,k) for k,v in cone_checks.items()}
height_sum=(I(delta)/pal_slope).hi
dx_bound=k*height_sum
dy_bound=delta+height_sum

# Given z in either overlap, its CB contact b has horizontal distance
# >= min_horizontal and total distance <=1, so a<= b_y-U_y <= vertical.
# Let r=b_y-C3_y.  The signed expansion is
#   a*r <= (b-U).(z-U)
#       <= A*k*d_y + vertical*delta - a*d_y.
# Averaging the two inequalities and using the height-sum bound gives
# a uniform bound on the mean r, without claiming that each r is small.
a=c3[1]-u[1]
assert a.lo>0 and (-m[5][1]).lo>0
horizontal=half-u[0]
min_horizontal=horizontal-dx_bound
assert min_horizontal.lo>0
vertical=(1-min_horizontal*min_horizontal).sqrt()+dy_bound
assert (horizontal*k-a).lo>0
r_bound=(( (horizontal*k-a)*I(height_sum/2) + vertical*delta )/a).hi
target=sub(c3,m[5])
normal,tangent=n[2],rot(n[2])
parallel=dot(normal,target)
transverse=dot(tangent,target)
lower_distance2=(parallel+half)*(parallel+half)+transverse*transverse
distance0=lower_distance2.sqrt()
midpoint_b=mul(-half,normal)
test_direction=mul(1/distance0,sub(target,midpoint_b))
direction_tangent=dot(test_direction,rot(m[5]))
tangent_coeff=max(abs(direction_tangent.lo),abs(direction_tangent.hi))
m5_slope=-m[5][1]
projected_error=((1+m5_slope)*I(r_bound)+I(tangent_coeff)*(2*I(r_bound)*m5_slope).sqrt()).hi
overlap_margin=distance0-(1+projected_error)
assert overlap_margin.lo>0,overlap_margin.out()

# Green integral over the additional cap C3 arc U1->U, old U->T chord
# and new T->U1 chord.  Its positive orientation is forced by convexity.
gain=(delta+cross(sub(c3,t),sub(u,u1)))/2
assert gain.lo>F(24,10**10), gain.out()
print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "coordinate_source_sha256":hashlib.sha256(raw).hexdigest(),
    "sigma_degrees_exact":str(sigma),
    "clockwise_U_shift_radians_exact":str(delta),
    "new_U_interval":[x.out() for x in u1],
    "geometry_strict_lower_bounds":geom_lower,
    "cone_strict_lower_bounds":cone_lower,
    "overlap_pair_total_height_below_new_U_upper":str(height_sum),
    "overlap_point_y_distance_from_old_U_upper":str(dy_bound),
    "overlap_points_x_displacement_above_original_U_upper":str(dx_bound),
    "mean_CB_contact_displacement_upper":str(r_bound),
    "mean_opposite_contact_projection_error_upper":str(projected_error),
    "projected_distance_overlap_contradiction_lower":str(overlap_margin.lo),
    "additional_cut_area_interval":gain.out(),
    "passed":True,
},indent=2))
