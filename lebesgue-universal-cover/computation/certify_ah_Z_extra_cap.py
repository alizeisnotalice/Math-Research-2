"""Quantitative new cap on the paired AH curve with certified continuous placement."""
from fractions import Fraction as F
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

source=Path(__file__).with_name("certify_ah_Z_root_rescue.py")
with contextlib.redirect_stdout(io.StringIO()) as stream:
    r=runpy.run_path(str(source))
prior=json.loads(stream.getvalue())
d=r["d"]
I,Jet,add,sub,mul,dot,cross,perp,cc=(d[k] for k in ("I","Jet","add","sub","mul","dot","cross","perp","cc"))
l0,mid,z,D=(d[k] for k in ("l0","mid","z","direction"))
half,r3=d["half"],d["r3"]
nl,nn,el,en=(tuple(map(I,d[k])) for k in ("n_l","n_n","e_l","e_n"))
trig=d["trig"]
old_t=d["t"]
epsilon=F(1,10**8)
delta=F(175,10**10)
w_cut=F(25,10**8)
co,si=(I(*trig.trig((delta,delta),k)) for k in ("cos","sin"))
d0=mul(-1,D)
length=dot(d0,d0).sqrt()
direct=mul(1/length,d0)
normal=perp(direct)
cut_dir=add(mul(co,direct),mul(si,normal))

def jet_vec(x):
    return tuple(map(Jet,x))

def paired_jets(w):
    jw=Jet(w)
    jw.c[1]=I(1)
    return d["pair"](add(jet_vec(z),mul(jw,jet_vec(D))))

def paired_point(w):
    return d["pair"](add(z,mul(I(w),D)))

q=paired_jets(I(0,w_cut))
vel=tuple(x.c[1] for x in q)
acc=tuple(2*x.c[2] for x in q)
speed=dot(vel,vel).sqrt()
assert speed.hi*w_cut<epsilon
assert dot(normal,acc).lo>0
# At zero the geometric pairing derivative is -Id.  The new chord
# enters the convex body along the retained side of the old straight.
assert si.lo>0 and co.lo>0
far=cross(cut_dir,sub(paired_point(w_cut),z))
assert far.lo>0
triangle=[paired_point(w_cut/4),paired_point(w_cut/2)]
removed=[-cross(cut_dir,sub(p,z)) for p in triangle]
assert all(x.lo>0 for x in removed)
gain=cross(sub(triangle[0],z),sub(triangle[1],z))/2
assert gain.lo>F(3,10**26)

# The inverse coordinate maps are roots of a unit-circle equation.
# Bound their gradients before using the mean-value theorem.
point_box=tuple(I(x.lo-epsilon,x.hi+epsilon) for x in z)
assert point_box[1].lo>F(518,1000)
assert dot(sub(point_box,d["c"]),sub(point_box,d["c"])).lo>1
def inverse_gradient(x,n,e):
    value=dot(x,n)-half
    root=(1-value*value).sqrt()
    return add(e,mul(value/root,n)),1/root

gr_l,lips_l=inverse_gradient(point_box,nl,el)
gr_n,lips_n=inverse_gradient(point_box,nn,en)
assert lips_l.hi<11 and lips_n.hi<4
# Verified two line roots remain on their small-root branches. The
# positivity of their square roots gives the orientation directly.
assert (point_box[1]-d["c"][1]).lo>0
# L,N chart center rectangle surrounds (ell=0,t=0) by more than 1e-5.
assert dot(sub(l0,d["c"]),el).lo>F(1,10**5)
assert dot(sub(d["g"],l0),el).lo>F(1,10**5)
assert dot(sub(d["e3"],mid),en).lo>F(1,10**5)
assert 11*epsilon<F(1,10**5)

# The positive old straight near Z supplies the graph g(ell). Its
# reflected branch follows by the exact pairing psi and t -> -t.
old_interval=add(z,mul(I(0,F(1,1000)),D))
old_l,_=inverse_gradient(old_interval,nl,el)
old_n,_=inverse_gradient(old_interval,nn,en)
ell_slope=dot(old_l,D)
nu_slope=dot(old_n,D)
assert ell_slope.hi<-F(11,1000) and nu_slope.lo>0
g_prime=nu_slope/ell_slope
assert -4<g_prime.lo<=g_prime.hi<-F(7,2)
ell_center=l0[0]+half+(r3/2)*I(F(-1,10**4),F(1,10**4))
psi_slope=half-r3/2*ell_center/(1-ell_center*ell_center).sqrt()
assert -F(21,20)<psi_slope.lo<=psi_slope.hi<-F(19,20)
assert F(19,20)*F(11,10**6)>F(1,10**5)
# Consequently -5 < g' < -3 on both branches of |ell|<=1e-5.

chart_l=I(F(-1,10**5),F(1,10**5))
chart_t=I(F(-1,10**5),F(1,10**5))
chart_lips=[]
for index in (0,1):
    a,b=Jet(chart_l),Jet(chart_t)
    (a if index==0 else b).c[1]=I(1)
    image=cc(add(jet_vec(l0),mul(a,jet_vec(el))),
             add(jet_vec(mid),mul(b,jet_vec(en))),-1)
    att=tuple(point.c[1] for point in image)
    chart_lips.append(dot(att,att).sqrt())
assert chart_lips[0].hi<F(1,5) and chart_lips[1].hi<F(3,5)

# Solid angle coefficients for the limiting pentagon. Each normalized
# tangent direction is a sum of the two endpoint normals with coefficient
# sum at most sec(half-angle) < 11/10.
sectors=r["order"]
for i,name in enumerate(sectors):
    a,b=r["edges"][name],r["edges"][sectors[(i+1)%len(sectors)]]
    assert dot(a,b).lo>F(3,4)
contact_error=60*epsilon
support_error=F(11,10)*contact_error
F_error=3*support_error
left=I(*map(F,prior["F_left_interval"]))
right=I(*map(F,prior["F_right_interval"]))
assert left.hi+F_error<0<right.lo-F_error
margin=F(prior["minimum_nonhex_edge_margin_lower"])
assert margin>5*support_error+2*delta
# R86 asserts margin > 0.000000821377744 at every nonhexagon edge.
# In particular its straight Ta-Z constraint tolerates this cap.
assert 2*delta<F("0.000000821377744")

print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "ideal_limit_source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
    "sigma_degrees_exact":prior["sigma_degrees_exact"],
    "baseline_EH_depth_radians_exact":prior["baseline_EH_depth_radians_exact"],
    "new_chord_angle_ccw_radians":str(delta),
    "paired_parameter_far_side_test":str(w_cut),
    "cap_radius_about_Z_upper":str(epsilon),
    "cap_triangle_area_lower_exact":str(gain.lo),
    "cap_triangle_area_interval":gain.out(),
    "paired_far_side_slack_lower":str(far.lo),
    "new_hat_side_triangle_slack_lower":[str(x.lo) for x in removed],
    "inverse_lipschitz_ell_upper":str(lips_l.hi),
    "inverse_lipschitz_t_upper":str(lips_n.hi),
    "original_graph_g_ell_derivative":g_prime.out(),
    "ell_pairing_derivative":psi_slope.out(),
    "chart_direct_derivative_norm_upper":[str(x.hi) for x in chart_lips],
    "double_hit_contact_error_upper_exact":str(contact_error),
    "support_error_upper_exact":str(support_error),
    "root_F_left_perturbed_upper":str(left.hi+F_error),
    "root_F_right_perturbed_lower":str(right.lo-F_error),
    "compatible_placement_clearance_after_cut_lower":str(margin-5*support_error-2*delta),
    "passed":True
},indent=2))
