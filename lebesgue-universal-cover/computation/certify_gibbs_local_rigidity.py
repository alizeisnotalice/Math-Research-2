"""Exact tests for the Reuleaux-triangle pinning and AH area density."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

source = Path(__file__).with_name("certify_gibbs_area_refined.py")
raw = source.read_bytes()
# The area source runs a whole quadrature after this boundary.  Here reuse only
# its reviewed interval, coordinate, and Taylor arithmetic prior to that loop.
parts = raw.decode().split("width = F(1,cells)")
assert len(parts) == 2
sys.argv = [str(source), "1.5494009336023418793", "2", "1"]
d = {"__file__": str(source)}
exec(compile(parts[0], str(source), "exec"), d)
I, Jet, add, sub, mul, dot, cross, cc = (
    d[key] for key in ("I","Jet","add","sub","mul","dot","cross","cc")
)
half, r3 = d["half"],d["r3"]
s3 = Jet(r3)
mid, z, t, n_l, n_n, e_l, e_n = (tuple(map(Jet,d[k])) for k in
                                  ("mid","z","t","n_l","n_n","e_l","e_n"))

# The exact wrong Reuleaux vertex is at polar angle 330 degrees and radius
# 1/sqrt(3).  The positive C cap normal is at 330+sigma degrees.
alpha = d["trig"].scale(d["pi_machin"], d["sigma"]/180)
bad_C_margin = I(*d["trig"].trig(alpha,"cos")) / r3 - half
assert bad_C_margin.lo > 0

# Exact convexification obstruction for the original EH implication (R33).
pal = (-r3/2, I(-half))
negative_e_arg = d["trig"].scale(d["pi_machin"], F(7,6)-d["sigma"]/180)
e_neg = tuple(I(*d["trig"].trig(negative_e_arg, kind)) for kind in ("cos","sin"))
h,u,te,ip,ss,b3,c3,e3,d1,f1,g,k = (d[key] for key in
    ("h","u","t_eh","ip","s","b3","c3","e3","d1","f1","g","k"))
rot = d["perp"]
eh_checks = {
    "EF_pal_F1_slack": half-dot(pal,f1),
    "EF_pal_H_slack": half-dot(pal,h),
    "HU_pal_H_derivative": dot(pal,rot(sub(h,c3))),
    "HU_pal_U_derivative": dot(pal,rot(sub(u,c3))),
    "TeI_negative_S_disk_T_derivative": -2*dot(sub(b3,ss),rot(sub(te,b3))),
    "TeI_negative_S_disk_I_derivative": -2*dot(sub(b3,ss),rot(sub(ip,b3))),
    "I_inside_S_disk": 1-dot(sub(ip,ss),sub(ip,ss)),
    "E3_inside_S_disk": 1-dot(sub(e3,ss),sub(e3,ss)),
    "E3_outside_negative_E_cap": half-dot(e_neg,e3),
    "D1_outside_negative_E_cap": half-dot(e_neg,d1),
    "right_side_E_cap_y_bound": half-F(3,5)/r3,
    "negative_E_cap_normal_x_negative": -e_neg[0],
    "negative_E_cap_normal_y_not_below_minus_3_5": F(3,5)+e_neg[1],
    "right_side_G_x_positive": g[0],
    "right_side_K_x_positive": k[0],
}
assert all(v.lo>0 for v in eh_checks.values())

def inverse(x, normal, along):
    dist = dot(normal,x)-half
    return sub(sub(x,mul(dist,normal)),mul((1-dist*dist).sqrt(),along))

def jacobian(x, l, n):
    w_l, w_n = sub(x,l),sub(x,n)
    lx = l[0]+half
    ds = (1-lx*lx).sqrt()*2/s3
    directional_l = dot(w_l,e_l)
    directional_n = dot(w_n,e_n)
    determinant = cross(w_l,w_n)
    # All factors except determinant are positive on the interior.  Only the
    # known corner P of the full rectangular map may have derivative zero.
    assert directional_l.c[0].lo>0 and directional_n.c[0].lo>0
    assert determinant.c[0].hi<0
    return -ds*directional_l*directional_n/determinant

def difference(position):
    w = Jet(position)
    w.c[1] = I(1)
    x = add(z,mul(w,sub(t,z)))
    l,n = inverse(x,n_l,e_l),inverse(x,n_n,e_n)
    c = l[0]+half
    lx = (c+s3*(1-c*c).sqrt())/2-half
    lp,np = (lx,(lx-1)/s3),sub(mul(2,mid),n)
    return jacobian(x,l,n)-jacobian(cc(lp,np,-1),lp,np)

origin = difference(I(0)).c[0]
assert origin.lo<=0<=origin.hi  # Equality algebraic: reflected parameters (0,0).
stack=[(F(0),F(1),0)]
leaves=[]
while stack:
    lo,hi,depth=stack.pop()
    try:
        slope=difference(I(lo,hi)).c[1]
        accepted=slope.lo>0
    except AssertionError:
        accepted=False
    if accepted:
        leaves.append((lo,hi,slope.lo))
        continue
    assert depth<18, (lo,hi)
    mid_w=(lo+hi)/2
    stack.extend(((lo,mid_w,depth+1),(mid_w,hi,depth+1)))
leaves.sort()
assert leaves[0][0]==0 and leaves[-1][1]==1
assert all(a[1]==b[0] for a,b in zip(leaves,leaves[1:]))
print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "area_source_sha256":hashlib.sha256(raw).hexdigest(),
    "sigma_degrees_exact":str(d["sigma"]),
    "excluded_C1_positive_depth_lower":str(bad_C_margin.lo),
    "AH_density_difference_derivative_minimum_lower":str(min(item[2] for item in leaves)),
    "AH_density_difference_derivative_partition":[[str(v) for v in leaf] for leaf in leaves],
    "AH_density_difference_at_Z_analytic":0,
    "EH_trigger_convexification_lower_bounds":{k:str(v.lo) for k,v in eh_checks.items()},
    "passed":True,
},indent=2))
