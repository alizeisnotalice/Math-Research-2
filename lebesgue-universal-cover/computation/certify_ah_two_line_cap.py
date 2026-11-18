"""Certify an AH envelope of two cuts at the deeper EH depth 0.00055."""
from fractions import Fraction as F
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

from flint import __version__ as flint_version, arb, ctx

AH_PRECISION = 128
ctx.prec = AH_PRECISION
here = Path(__file__).resolve().parent
argv = sys.argv[:]
sigma_degrees = "1.5497"
eh_depth = "0.00055"
with contextlib.redirect_stdout(io.StringIO()):
    audit = runpy.run_path(str(here / "ah_root_geometry_at_angle.py"),
                           init_globals={"sigma_degrees": sigma_degrees, "eh_depth": eh_depth})
sys.argv = argv
d = audit["d"]


def A(q):
    if isinstance(q, arb):
        return q
    if isinstance(q, d["I"]):
        return box(q.lo, q.hi)
    q = F(q)
    return arb(q.numerator) / q.denominator


def box(lo, hi):
    a,b=A(lo),A(hi)
    return a.union(b)


def av(p):
    return tuple(map(A,p))


add,sub,mul,dot,cross,perp = (d[k] for k in ("add","sub","mul","dot","cross","perp"))
z,D,l0,mid = (av(d[k]) for k in ("z","direction","l0","mid"))
nl,nn,el,en = (av(d[k]) for k in ("n_l","n_n","e_l","e_n"))
r3,half = A(d["r3"]),A("0.5")
def norm(v):
    return dot(v,v).sqrt()
length = norm(D)
direct = mul(-1/length,D)
normal = perp(direct)
zero,one = A(0),A(1)
cuts=[(F("0.355"), F("0.00060")), (F("0.410"), F("0.00030"))]
wp, wf, wmax, eps = F("0.410"), F("0.07"), F("0.43"), F("0.00065")
pbound = F("0.072")
tangents=[A(t).tan() for _,t in cuts]
pivot_set=[add(z,mul(A(p),D)) for p,_ in cuts]
cut_dirs=[add(mul(A(a).cos(),direct),mul(A(a).sin(),normal)) for _,a in cuts]

def height_bound(w):
    return max((t*length*A(p-w)).upper() for (p,_),t in zip(cuts,tangents))


def inverse(x,n,e):
    h=dot(n,x)-half
    root=(1-h*h).sqrt()
    return sub(sub(x,mul(h,n)),mul(root,e)), add(e,mul(h/root,n))

def cc(a,b):
    diff=sub(b,a)
    return sub(mul(half,add(a,b)),mul((1/dot(diff,diff)-A("0.25")).sqrt(),perp(diff)))

def X(ell,t):
    return cc(add(l0,mul(ell,el)),add(mid,mul(t,en)))

def derivs(x,l,n):
    a,b=sub(x,l),sub(x,n)
    det=cross(a,b)
    xl=mul(-dot(a,el)/det,perp(b))
    xt=mul(dot(b,en)/det,perp(a))
    return xl,xt

def psi_prime(l):
    u=l[0]+half
    return half-r3*half*u/(1-u*u).sqrt()

def reflection(v):
    return (half*(r3*v[1]-v[0]),half*(r3*v[0]+v[1]))

def ell_partner(l):
    lx = l[0]+half
    xx=(lx+r3*(1-lx*lx).sqrt())/2-half
    return (xx,(xx-1)/r3)

def points(w):
    x=add(z,mul(w,D))
    l,grad_l=inverse(x,nl,el)
    n,grad_n=inverse(x,nn,en)
    lp=ell_partner(l)
    return (l,reflection(lp),reflection(cc(lp,sub(mul(2,mid),n))),n,x),grad_l,grad_n

def point_jets(w):
    vv,gl,gn=points(w)
    l,lp,xp,n,x=vv
    lp,xp=reflection(lp),reflection(xp)
    np=sub(mul(2,mid),n)
    dl,dn=dot(gl,D),dot(gn,D)
    dlp=psi_prime(l)*dl
    xp_l,xp_t=derivs(xp,lp,np)
    xp_w=sub(mul(dlp,xp_l),mul(dn,xp_t))
    return vv, (mul(dl,el),reflection(mul(dlp,el)),reflection(xp_w),mul(dn,en),D)

def centered_points(a,b):
    middle=points(A((a+b)/2))[0]
    vv,prime=point_jets(box(a,b))
    step=box((a-b)/2,(b-a)/2)
    return tuple(add(val,mul(step,slope)) for val,slope in zip(middle,prime))

def q(w):
    return reflection(points(A(w))[0][2])

for pivot,c in zip(pivot_set,cut_dirs):
    assert cross(c,sub(q(wf),pivot))>0
    assert cross(c,sub(z,pivot))<0

# R86 verifies these coordinates and the retained boundary decomposition.
# Recheck the endpoint signs for the five exact Sprague circular pieces.
hp = list(map(av,d["n"]))
mp = list(map(av,d["m"]))
coords = {k: av(d[k]) for k in
          ("b1","b3","c3","d1","e3","f1","f3","g","h","ip","k","p","q","t_eh")}
def meet(a,b):
    delta=cross(a,b)
    return (half*(b[1]-a[1])/delta,half*(a[0]-b[0])/delta)
def get(name):
    return coords[name]
ZP=reflection(points(one)[0][2])
b2,f2,c2 = (meet(hp[j],mp[k]) for j,k in ((0,0),(2,2),(5,5)))
e1=meet(hp[3],hp[4])
arc_tests=[
    ("Q-Y",get("e3"),get("d1"),get("q"),av(d["y"])),
    ("Tprime-P",get("g"),get("d1"),get("p"),ZP),
    ("H-U1",get("c3"),get("b1"),get("h"),av(audit["u1"])),
    ("Te-I",get("b3"),b2,get("ip"),get("t_eh")),
    ("G-K",get("f3"),f2,get("k"),get("g")),
]
arc_signs={}
for label,center,other,tangent,strict in arc_tests:
    tangent_radial,far_radial=sub(tangent,center),sub(strict,center)
    chord=sub(other,center)
    # The tangent endpoint equality follows from orthogonal projection.
    assert dot(tangent_radial,chord).contains(0)
    assert dot(tangent_radial,far_radial)>0
    sign=-dot(far_radial,chord)
    assert sign>0, (label,sign)
    arc_signs[label]=str(sign.lower())
for start,end,n in ((get("b1"),get("c3"),mp[2]),(e1,get("f1"),mp[0])):
    assert dot(start,n)<half and dot(end,n)<half
# This is the DC endpoint implication that yields [D1,G].
v=sub(c2,get("f3"))
ug=sub(get("g"),get("f3"))
assert dot(ug,v)>1 and dot(v,v)>1
assert dot(v,sub(f2,get("f3")))<0


def chart(x):
    l,_=inverse(x,nl,el)
    n,_=inverse(x,nn,en)
    return dot(sub(l,l0),el),dot(sub(n,mid),en)

def psi(ell):
    l=add(l0,mul(ell,el))
    return dot(sub(ell_partner(l),l0),el)

def straight_graph(ell):
    l=add(l0,mul(ell,el))
    radius=dot(normal,sub(z,l))
    # The upper normal branch has a positive projection on direct.
    x=add(add(l,mul(radius,normal)),mul((1-radius*radius).sqrt(),direct))
    return chart(x)[1]

def graph_value(l):
    return straight_graph(l) if l<=0 else -straight_graph(psi(l))

def graph_hull(l):
    return coord_hull(graph_value,l)

def chart_bounds(ell,t):
    x=X(ell,t)
    l=add(l0,mul(ell,el))
    n=add(mid,mul(t,en))
    return x,l,n

def coord_hull(f,ell):
    # Both call sites are monotone on the specified interval.
    return f(ell.lower()).union(f(ell.upper()))

assert dot(direct,sub(z,l0))>0
assert abs(straight_graph(0))<A(F(1,10**30))

for j in range(int(1000*wf)):
    qa=reflection(centered_points(F(j,1000),F(j+1,1000))[2])
    assert dot(direct,sub(qa,z))/length<A(pbound)
hmax=height_bound(-pbound)
assert hmax<A("0.000012")
lmin,tmax=chart(add(z,mul(A(wmax),D)))
lmax=psi(lmin)
max_xl,max_xt=zero,zero
min_eta_speed=one
psi_max=zero
for j in range(int(1000*wmax)):
    a,b=F(j,1000),F(j+1,1000)
    wa,wb=chart(add(z,mul(A(a),D))),chart(add(z,mul(A(b),D)))
    l=wa[0].union(wb[0])
    t=wa[1].union(wb[1])
    _v,gl,gn=points(box(a,b))
    assert dot(gl,D)<0<dot(gn,D)
    slope=dot(gn,D)/dot(gl,D)
    assert slope<A("-3.5")
    p_l=coord_hull(psi,l)
    for E,T in ((l,t),(p_l,-t)):
        ps=-psi_prime(add(l0,mul(E,el)))
        assert ps>A("0.95") and ps<A("1.05")
        psi_max=max(psi_max,ps.upper())
        ln=E + box(-eps,eps)
        tn=T + box(-eps,eps)
        x,lc,nc=chart_bounds(ln,tn)
        xl,xt=derivs(x,lc,nc)
        max_xl,max_xt=max(max_xl,norm(xl).upper()),max(max_xt,norm(xt).upper())
        norm_l, norm_t = ("0.24", "0.68")
        assert norm(xl) < A(norm_l) and norm(xt) < A(norm_t), (a, norm(xl), norm(xt))

gate_C=av(d["c"])
left_l=dot(sub(gate_C,l0),el)
right_l=dot(sub(get("g"),l0),el)
tau=dot(sub(get("e3"),mid),en)
ref_l=[chart(add(z,mul(A(F(i,100)),D)))[0] for i in range(int(100*wmax)+1)]
ref_boxes=[a.union(b) for a,b in zip(ref_l,ref_l[1:])]
deficit_by_bin = [zero for _ in ref_boxes]
parameter_boxes=[]
for j in range(int(1000*(wp+pbound))):
    a,b=-pbound+F(j,1000),-pbound+F(j+1,1000)
    w=box(a,b)
    own_height = height_bound(a)
    x=add(add(z,mul(w,D)),mul(box(0,own_height.upper()),normal))
    cl,ct=chart(add(z,mul(A((a+b)/2),D)))
    _,gl=inverse(x,nl,el)
    _,gt=inverse(x,nn,en)
    step=box((a-b)/2,(b-a)/2)
    ht=box(0,own_height.upper())
    l=cl+dot(gl,D)*step+dot(gl,normal)*ht
    t=ct+dot(gt,D)*step+dot(gt,normal)*ht
    assert left_l<l<right_l and -tau<t<tau
    lc=add(l0,mul(l,el))
    nc=add(mid,mul(t,en))
    assert cross(sub(nc,lc),sub(x,lc))<0

    pl=coord_hull(psi,l)
    assert x[1]>A("0.518") and dot(sub(x,gate_C),sub(x,gate_C))>1,(a,x)
    g=graph_hull(l)
    xx,lc,nc=chart_bounds(l,t.union(g))
    xl,xt=derivs(xx,lc,nc)
    value=-dot(normal,xt)
    min_eta_speed = min(min_eta_speed, value.lower())
    assert value > A("0.008"), (a, value)
    cell_deficit = (own_height / value.lower()).upper()
    assert cell_deficit < A(eps)
    # Pairing sends whichever of l and psi(l) is nonpositive to the
    # straight graph. If the sign is unresolved, test both interval images.
    used = {k for v in (l, pl) if not v > 0
            for k, ell_range in enumerate(ref_boxes) if v.overlaps(ell_range)}
    assert used
    parameter_boxes.append((l,pl,used,a))
    for k in used:
        deficit_by_bin[k] = max(deficit_by_bin[k], cell_deficit)

for l,pl,used,a in parameter_boxes:
    margin=A(2)/3*max(deficit_by_bin[k] for k in used)
    assert l>lmin+margin and l<lmax-margin,(a,l,margin)
    assert pl>lmin+margin and pl<lmax-margin,(a,pl,margin)



# A bracket reference graph with round-to-0.001 degree centres. The
# reference is merely a choice of angles, checked below by Arb.
root_seeds=[19844,19565,19296,19035,18783,18537,18298,18065,17837,17612,
            17390,17163,16935,16703,16465,16224,15978,15732,15481,15229,
            14974,14713,14448,14178,13905,13626,13343,13054,
            12760,12460,12154,11840,11519,11189,10850,10500,10136,9756,9357,8932,8469,7948]
low_root_seeds=[None]*31+[3326,3464,3604,3746,3890,4036,4183,4332,4484,4657,4886,5231]
orders=audit["order"]
sector=audit["sector_formula"]
rows=[(label,av(v),A(c)) for label,v,c in audit["rows"]
      if label.startswith(("Y-Ta","Ta-Z","U1-Te","I-E3","K-C3"))]
previous=z
for i in range(1,17):
    after=q(F(i,16))
    normal_to_chord=mul(-1/norm(sub(after,previous)),perp(sub(after,previous)))
    rows.append((f"Z-Tprime:{i}",normal_to_chord,dot(normal_to_chord,previous)))
    previous=after
for j,(cut_direction,pivot) in enumerate(zip(cut_dirs,pivot_set)):
    cap_normal=mul(-1,perp(cut_direction))
    rows.append((f"AH-line:{j}",cap_normal,dot(cap_normal,pivot)))
n0,n1=map(av,d["n"][:2])
angle_unit = arb.pi()/180
def trig(phi):
    a=phi*angle_unit
    return (a.cos(),a.sin())
def rotate(v,sincos):
    return add(mul(sincos[0],v),mul(sincos[1],perp(v)))
def supports(v,e,sincos,pre,differentiated=False):
    u=rotate(pre,sincos)
    pp=[]
    cr=[cross(x,u) for x in e]
    for j,(k,c) in enumerate(sector):
        if not cr[j]<0 and not cr[(j+1)%10]>0:
            pp.append(dot(v[k],perp(u)) if differentiated else dot(v[k],u)+c)
    assert pp, cr
    ans=pp[0]
    for p in pp[1:]:
        ans=ans.union(p)
    return ans

def F_at(v,e,phi):
    ans=-half
    for j in range(3):
        ans+=(-1)**j*supports(v,e,trig(A(phi+60*j)),n0)
    return ans

def placement(v,e,a,b,turn,sign):
    mt=trig(A((a+b)/2-turn))
    sn=trig(box(a-turn,b-turn))
    dw=angle_unit*box((a-b)/2,(b-a)/2)
    n1p=(n1[0],sign*n1[1])
    a0=half-supports(v,e,mt,n0)
    a1=half-supports(v,e,mt,n1p)
    da0=-supports(v,e,sn,n0,True)
    da1=-supports(v,e,sn,n1p,True)
    trans=(a0,(2*a1-a0)/r3)
    dtrans=(da0,(2*da1-da0)/r3)
    margin=[]
    for _,u,c in rows:
        u_ref=(u[0],sign*u[1])
        at_mid=c-supports(v,e,mt,u_ref)-dot(u,trans)
        derivative=-supports(v,e,sn,u_ref,True)-dot(u,dtrans)
        margin.append(at_mid+dw*derivative)
    return margin

angular_cells=32
counts=0
max_gap=zero
best_root=arb(1)
best_edge=arb(1)
best_actual_edge=arb(1)
for j,local_error in enumerate(deficit_by_bin):
    if local_error==zero:
        continue
    if j<9:
        width=("2.4","2.4","2.4","2.4","2.4","2.3","2.1","1.9","1.75")[j]
        turn,sign,half_bracket,cells = 300,1,F(width),40
    elif j<30:
        width="1.7" if j<26 else ("1.6","1.6","1.4","1.2")[j-26]
        turn,sign,half_bracket,cells = 120,-1,F(width),40
    elif j<35:
        width=("1","1","1","0.9","0.8")[j-30]
        turn,sign,half_bracket,cells = 300,1,F(width),75
    else:
        width=("1.9","2","2.1","2.15","1.65","1.25","0.65")[j-35]
        turn,sign,half_bracket,cells = 0,1,F(width),90
    local_seeds=root_seeds if turn else low_root_seeds
    # All five points have distance at most the local deficit. The DE
    # contact is on a unit-coordinate line and can attain that full error.
    support_error = A(F(107,100))*local_error
    gap = 3*support_error
    max_gap = max(max_gap,gap.upper())

    for i in range(cells):
        counts+=1
        a=F(j,100)+F(i,100*cells)
        b=F(j,100)+F(i+1,100*cells)
        w=box(a,b)
        vertices=centered_points(a,b)
        _,grad_l,grad_n=points(w)
        edges=[sub(vertices[k],vertices[l]) for k,l in orders]
        for k,x in enumerate(edges):
            y=edges[(k+1)%10]
            assert (x[1]>0 if k<5 else x[1]<0), (a,k,x)
            assert cross(x,y)>0 and dot(x,y)>A("0.75"), (a,k,x,y)
        for k in range(5):
            for l in range(k+1,5):
                if (k,l) not in orders and (l,k) not in orders:
                    e=sub(vertices[k],vertices[l])
                    assert dot(e,e)<1
        assert dot(grad_n,D)/dot(grad_l,D)<A("-3.5"), (a,grad_l,grad_n)
        wcen=F(2*i+1,2*cells)
        predict=((1-wcen)*local_seeds[j]+wcen*local_seeds[j+1])/1000
        al,bl=predict-half_bracket,predict+half_bracket
        sgn=1 if turn else -1
        fs=[sgn*F_at(vertices,edges,al)-gap,
            -sgn*F_at(vertices,edges,bl)-gap]
        assert min(f.lower() for f in fs)>0, (a,al,bl,fs)
        best_root=min(best_root,*(f.lower() for f in fs))
        for p in range(angular_cells):
            angle_left=al+(bl-al)*F(p,angular_cells)
            angle_right=al+(bl-al)*F(p+1,angular_cells)
            ms=placement(vertices,edges,angle_left,angle_right,turn,sign)
            worst=min(range(len(ms)),key=lambda k:ms[k].lower())
            assert ms[worst]>A(gap), (a,p,rows[worst][0],ms[worst].lower(),gap)
            best_edge=min(best_edge,ms[worst].lower())
            best_actual_edge=min(best_actual_edge,(ms[worst]-gap).lower())


eh_path = here / "eh-cap-switch-joint-550.json"
baseline = json.loads(eh_path.read_bytes())["summary"]
assert F(baseline["sigma_degrees"]) == F(sigma_degrees)
assert F(baseline["delta"]) == F(eh_depth)
assert baseline["validated_duals"] == 60 and baseline["EH_combined_cells"] > 0

print(json.dumps({
    "python":sys.version,
    "python_flint":flint_version,
    "arb_precision_bits":AH_PRECISION,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "geometry_sha256":hashlib.sha256((here/"ah_root_geometry_at_angle.py").read_bytes()).hexdigest(),
    "EH_baseline_sha256":hashlib.sha256(eh_path.read_bytes()).hexdigest(),
    "sigma_degrees_exact":sigma_degrees,
    "EH_depth_radians_exact":eh_depth,
    "AH_cut_parameters_exact":[[str(p),str(a)] for p,a in cuts],
    "AH_paired_arc_far_parameter_exact":str(wf),
    "AH_cap_deficit_by_reference_bin":[str(value) for value in deficit_by_bin],
    "AH_support_parameter_cells":counts,
    "angle_cells_per_compatibility_bracket":angular_cells,
    "uniform_perturbed_compatibility_slack_lower":str(best_root),
    "uniform_actual_placement_slack_lower":str(best_actual_edge),
    "AH_normal_speed_lower":str(min_eta_speed),
    "AH_graph_partial_derivative_norm_uppers":[str(max_xl),str(max_xt)],
    "tested_support_rows":[r[0] for r in rows],
    "true_Sprague_endpoint_signs":arc_signs,
    "passed":True,
},indent=2))
