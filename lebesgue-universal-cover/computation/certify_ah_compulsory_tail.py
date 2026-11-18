"""Certify compulsory AH tails and the true G arc for a continuum of parameters."""
from fractions import Fraction as F
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

from flint import arb, ctx

ctx.prec = 160
here = Path(__file__).resolve().parent
source = here / "ah_root_geometry_at_angle.py"
with contextlib.redirect_stdout(io.StringIO()):
    audit = runpy.run_path(str(source),
                          init_globals={"sigma_degrees":"1.5497","eh_depth":"0.00036"})
d = audit["d"]
I = d["I"]
add, sub, mul, dot, cross, perp = (d[key] for key in ("add", "sub", "mul", "dot", "cross", "perp"))


def ball(x):
    if isinstance(x, arb):
        return x
    if isinstance(x, I):
        return ball(x.lo).union(ball(x.hi))
    x = F(x)
    return arb(x.numerator)/x.denominator


def box(a,b):
    return ball(a).union(ball(b))


def av(v):
    return tuple(map(ball,v))


z, D, l0, M, n_l, e_l, n_n, e_n = (av(d[k]) for k in (
    "z", "direction", "l0", "mid", "n_l", "e_l", "n_n", "e_n"))
r3, half = ball(d["r3"]), ball(F(1,2))
n, m = (list(map(av,d[k])) for k in ("n","m"))
C = av(d["c"])
T = av(d["t"])
m5_prime = (-m[3][0],m[3][1])
assert dot(m5_prime,e_l) > 0
for x in (add(z,mul(ball("0.429"),D)),T):
    assert dot(sub(x,C),e_l) > 0


def inverse(x, normal, along):
    offset = dot(x,normal)-half
    root = (1-offset*offset).sqrt()
    center = sub(sub(x,mul(offset,normal)),mul(root,along))
    grad = add(along,mul(offset/root,normal))
    return center, grad


def reflect(v):
    return (half*(r3*v[1]-v[0]), half*(r3*v[0]+v[1]))


def intersection(a,b):
    q = sub(b,a)
    return sub(mul(half,add(a,b)),mul((1/dot(q,q)-ball(F(1,4))).sqrt(),perp(q)))


def direct_derivatives(x,l,N):
    a,b = sub(x,l),sub(x,N)
    determinant = cross(a,b)
    return mul(-dot(a,e_l)/determinant,perp(b)),mul(dot(b,e_n)/determinant,perp(a))


def pattern(w):
    x = add(z,mul(w,D))
    L, grad_l = inverse(x,n_l,e_l)
    N, grad_n = inverse(x,n_n,e_n)
    a = L[0]+half
    radius = (1-a*a).sqrt()
    lx_partner = (a+r3*radius)/2-half
    Lp = (lx_partner,(lx_partner-1)/r3)
    Np = sub(mul(2,M),N)
    xp = intersection(Lp,Np)
    dl,dn = dot(grad_l,D),dot(grad_n,D)
    dp = (half-r3*half*a/radius)*dl
    px,pt = direct_derivatives(xp,Lp,Np)
    dxp = sub(mul(dp,px),mul(dn,pt))
    return ((L,reflect(Lp),reflect(xp),N,x),
            (mul(dl,e_l),reflect(mul(dp,e_l)),reflect(dxp),mul(dn,e_n),D))


order, formulas = audit["order"], audit["sector_formula"]
for j,(k,c) in enumerate(formulas):
    for i in (j,(j+1)%10):
        assert order[i][1 if c else 0] == k


def sectors(vertices):
    return [sub(vertices[i],vertices[j]) for i,j in order]


def choices(rays,u):
    signs = [cross(r,u) for r in rays]
    return [(k,c) for i,(k,c) in enumerate(formulas)
            if not signs[i]<0 and not signs[(i+1)%10]>0]


def hull(values):
    it=iter(values)
    result=next(it)
    for x in it:
        result=result.union(x)
    return result


def support(vertices,rays,u, derivative=None):
    opt = choices(rays,u)
    if derivative=="angle":
        return hull(dot(vertices[k],perp(u)) for k,c in opt)
    if derivative=="parameter":
        return hull(dot(vertices[k],u) for k,c in opt)
    return hull(dot(vertices[k],u)+c for k,c in opt)


radian = arb.pi()/180


def unit(alpha):
    a=ball(alpha)*radian
    return a.cos(),a.sin()


def compatibility(vertices,rays,alpha,derivative=None):
    v = sum((-1)**j*support(vertices,rays,unit(alpha+60*j),derivative) for j in range(3))
    if derivative=="angle":
        return v*radian
    return v if derivative else v-half


def inverse_pose(u,turn,sgn):
    w=(u[0],sgn*u[1])
    co,si=unit(-turn)
    return add(mul(co,w),mul(si,perp(w)))


allowed = {(0,1),(240,-1)}
vertical = (180,-1)
other_poses = [(turn,sgn) for turn in range(0,360,60) for sgn in (1,-1)
               if (turn,sgn) not in allowed | {vertical}]
cap_tests=[[inverse_pose(mj,turn,sgn) for mj in (m[3],m[5])]
           for turn,sgn in other_poses]

G, E3 = av(d["g"]), av(d["e3"])
N_at_T,_ = inverse(T,n_n,e_n)
t_begin = dot(sub(N_at_T,M),e_n)
t_end = dot(sub(E3,M),e_n)
t_velocity = t_end - t_begin
assert t_velocity > 0 and dot(m[3],e_n)>0

def circle_pattern(parameter):
    t = t_begin + parameter*t_velocity
    N = add(M,mul(t,e_n))
    N_other = sub(mul(2,M),N)
    x = intersection(C,N)
    x_other = intersection(G,N_other)
    _,xp = direct_derivatives(x,C,N)
    _,xp_other = direct_derivatives(x_other,G,N_other)
    zero=(ball(0),ball(0))
    return ((C,reflect(G),reflect(x_other),N,x),
            (zero,zero,mul(-t_velocity,reflect(xp_other)),mul(t_velocity,e_n),mul(t_velocity,xp)))


border = F("0.05")

def verify_grid(start,end,width,geometry,caps,circle_arc=False):
    min_root,min_pose = ball(1),ball(1)
    left_slope,right_slope = -ball(1),ball(1)
    tested = successful_cells = 0
    a = start
    while a < end:
        b=min(a+width,end)
        wmid=(a+b)/2
        middle,_=geometry(ball(wmid))
        _,w_derivative=geometry(box(a,b))
        dw=box(a-wmid,b-wmid)
        vertices=[add(v,mul(dw,vprime)) for v,vprime in zip(middle,w_derivative)]
        rays=sectors(vertices)
        rays_mid=sectors(middle)
        for i,p in enumerate(rays):
            q=rays[(i+1)%10]
            if circle_arc and i in (0,5):
                # The exact DC lower-root sign pins these antipodal
                # rays to the appropriate closed halfplanes at P.
                assert (p[0]>0 if i==0 else p[0]<0), (a,i,p)
            else:
                assert (p[1]>0 if i<5 else p[1]<0), (a,i,p)
            assert cross(p,q)>0 and dot(p,q)>0,(a,i,p,q)
        for i in range(5):
            for j in range(i+1,5):
                if j-i not in (1,4):
                    delta=sub(vertices[i],vertices[j])
                    assert dot(delta,delta)<1,(a,i,j)
        for u,k in ((n[4],3),(n[5],0)):
            assert choices(rays,u)==[(k,0)],(a,k,choices(rays,u))
        if circle_arc:
            # On the G arc the last ray of the EF vertex reaches n3
            # exactly at P. The remaining sign is exact by the lower
            # DC root: cross(n3,v[1]-v[2])=e_l.(x_other-G)>=0.
            assert cross(rays[4],n[3]) > 0
            if a>=F("0.85"):
                assert choices(rays,m[3]) == [(3,0)], (a,choices(rays,m[3]))
            else:
                assert support(vertices,rays,m[0]) > half
        else:
            assert choices(rays,n[3]) == [(1,0)]
        for cases in caps:
            loss=max((support(vertices,rays,u)-half).lower() for u in cases)
            assert loss>0,(a,loss)
            min_pose=min(min_pose,loss)
        fl=compatibility(vertices,rays,box(0,border),"angle")
        fr=compatibility(vertices,rays,box(60-border,60),"angle")
        assert fl<0<fr,(a,fl,fr)
        left_slope=max(left_slope,fl.upper())
        right_slope=min(right_slope,fr.lower())

        def certify_angle(lo,hi,level=0):
            nonlocal tested,min_root
            cm=(lo+hi)/2
            phi=box(lo,hi)
            first=compatibility(middle,rays_mid,ball(cm))
            slope_angle=compatibility(vertices,rays,phi,"angle")
            slope_w=compatibility(w_derivative,rays,phi,"parameter")
            value=first+box(lo-cm,hi-cm)*slope_angle+dw*slope_w
            if value<0:
                tested+=1
                min_root=min(min_root,(-value).lower())
                return
            assert level<9,(a,b,lo,hi,value)
            certify_angle(lo,cm,level+1)
            certify_angle(cm,hi,level+1)

        delta_phi=F(1,4)
        lo=border
        while lo<60-border:
            hi=min(lo+delta_phi,60-border)
            certify_angle(lo,hi)
            lo=hi
        successful_cells+=1
        a=b
    return {
        "parameter_cells":successful_cells,
        "accepted_joint_cells":tested,
        "F_left_angle_slope_upper":str(left_slope),
        "F_right_angle_slope_lower":str(right_slope),
        "F_midrange_negative_slack_lower":str(min_root),
        "uniform_slanted_support_violation_lower":str(min_pose),
    }


w_start,w_end=F("0.429"),F(1)
tail = verify_grid(w_start,w_end,F("0.002"),pattern,cap_tests)
# The rotation R300 is excluded on the G-arc family by m3·N(t)<1/2
# and h(m0)=1-h(m3), with N(t) the exact m3 support.
arc_caps = [cases for key,cases in zip(other_poses,cap_tests) if key!=(300,1)]
arc = verify_grid(F(0),F(1),F("0.01"),circle_pattern,arc_caps,circle_arc=True)

length = dot(D,D).sqrt()
direct = mul(-1/length,D)
normal = perp(direct)
projection_speed = ball(1)
for j in range(1000):
    _,velocities = pattern(box(F(j,1000),F(j+1,1000)))
    speed = dot(direct,reflect(velocities[2]))
    assert speed > 0, (j,speed)
    projection_speed = min(projection_speed,speed.lower())
at_start,_ = pattern(ball(w_start))
from_z = sub(reflect(at_start[2]),z)
projected_lower = dot(direct,from_z)/length
height_lower = dot(normal,from_z)/length
assert projected_lower > 0 and height_lower > 0
g_arc_start = reflect(pattern(ball(1))[0][2])
P = av(d["p"])
assert dot(direct,perp(sub(g_arc_start,G))) > 0
assert dot(direct,perp(sub(P,G))) > 0
max_slope_p0 = height_lower/projected_lower
max_slope_p41 = height_lower/(ball("0.410")+projected_lower)

print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "geometry_source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
    "sigma_degrees_exact":"1.5497",
    "compulsory_tail_parameter_range_exact":[str(w_start),str(w_end)],
    "angle_endpoint_collar_degrees":str(border),
    "straight_and_paired_tail":tail,
    "true_G_arc":arc,
    "analytic_pose_exclusions":["vertical reflection by exact DC contact gate",
                                "R300 on G arc by exact DE support at m3"],
    "q_projection_speed_lower":str(projection_speed),
    "forced_arc_endpoint_normalized_forward_projection":str(projected_lower),
    "forced_arc_endpoint_normalized_height":str(height_lower),
    "necessary_max_slope_at_p0":str(max_slope_p0),
    "necessary_max_slope_at_p410":str(max_slope_p41),
    "passed":True,
},indent=2))
