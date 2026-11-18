"""Arb certificate for a finite AH cut and its limiting pentagon placement."""
from fractions import Fraction as F
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

from flint import __version__ as flint_version, arb, ctx

ctx.prec = 128
here = Path(__file__).resolve().parent
argv = sys.argv[:]
with contextlib.redirect_stdout(io.StringIO()):
    audit = runpy.run_path(str(here / "certify_ah_Z_root_rescue.py"))
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
delta = F("0.00025")
wp = F(1,4)
wf = F(1,25)
wmax = F(27,100)
eps = F("0.00021")
cut_d=add(mul(A(delta).cos(),direct),mul(A(delta).sin(),normal))
pivot=add(z,mul(A(wp),D))

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

def cut(x):
    return cross(cut_d,sub(x,pivot))

assert cut(q(wf))>0 and cut(z)<0
allowed = A(delta).tan()*length*A(wp)
height = allowed * A("0.99")
edge_q = q(wf)
assert dot(direct,sub(edge_q,z))>0
height_q = A(wp)*length*dot(normal,sub(edge_q,z))/dot(direct,sub(edge_q,pivot))
assert height_q>height
s=add(z,mul(height,normal))
assert cut(s)<0
gain=cross(sub(z,pivot),sub(s,pivot))/2
gain_lower = F("0.00000001158419305978677")
assert gain>A(gain_lower)

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

pbound=F(41,1000)
for j in range(40):
    qa=reflection(centered_points(F(j,1000),F(j+1,1000))[2])
    assert dot(direct,sub(qa,z))/length<A(pbound)
hmax=A(delta).tan()*length*A(wp+pbound)
assert hmax<A("0.000002816")
lmin,tmax=chart(add(z,mul(A(wmax),D)))
lmax=psi(lmin)
buffer=A(2*eps/3)
max_xl,max_xt=zero,zero
min_eta_speed=one
psi_max=zero
for j in range(270):
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
        assert norm(xl)<A("0.2") and norm(xt)<A("0.6")

gate_C=av(d["c"])
left_l=dot(sub(gate_C,l0),el)
right_l=dot(sub(get("g"),l0),el)
tau=dot(sub(get("e3"),mid),en)
for j in range(291):
    a,b=-pbound+F(j,1000),-pbound+F(j+1,1000)
    w=box(a,b)
    x=add(add(z,mul(w,D)),mul(box(0,hmax.upper()),normal))
    cl,ct=chart(add(z,mul(A((a+b)/2),D)))
    _,gl=inverse(x,nl,el)
    _,gt=inverse(x,nn,en)
    step=box((a-b)/2,(b-a)/2)
    ht=box(0,hmax.upper())
    l=cl+dot(gl,D)*step+dot(gl,normal)*ht
    t=ct+dot(gt,D)*step+dot(gt,normal)*ht
    assert l>lmin+buffer and l<lmax-buffer, (a,l,lmin,lmax)
    assert left_l<l<right_l and -tau<t<tau
    lc=add(l0,mul(l,el))
    nc=add(mid,mul(t,en))
    assert cross(sub(nc,lc),sub(x,lc))<0

    pl=coord_hull(psi,l)
    assert pl>lmin+buffer and pl<lmax-buffer,(a,pl,lmin,lmax)
    assert x[1]>A("0.518") and dot(sub(x,gate_C),sub(x,gate_C))>1,(a,x)
    g=graph_hull(l)
    xx,lc,nc=chart_bounds(l,t.union(g))
    xl,xt=derivs(xx,lc,nc)
    value=-dot(normal,xt)
    assert value>A("0.014"), (a,value)
    min_eta_speed=min(min_eta_speed,value.lower())
assert hmax<A(eps)*A("0.014")


# A bracket reference graph with round-to-0.001 degree centres. The
# reference is merely a choice of angles, checked below by Arb.
root_seeds=[19844,19565,19296,19035,18783,18537,18298,18065,17837,17612,
            17390,17163,16935,16703,16465,16224,15978,15732,15481,15229,
            14974,14713,14448,14178,13905,13626,13343,13054]
orders=audit["order"]
sector=audit["sector_formula"]
rows=[(label,av(v),A(c)) for label,v,c in audit["rows"]
      if label.startswith(("Y-Ta","Ta-Z","Z-Tprime","U1-Te","I-E3","K-C3"))]
new_cap_normal=mul(-1,perp(cut_d))
rows.append(("AH-cap",new_cap_normal,dot(new_cap_normal,pivot)))
n0,n1=map(av,d["n"][:2])
angle_unit = arb.pi()/180
def trig(phi):
    a=phi*angle_unit
    return (a.cos(),a.sin())
def rotate(v,sincos):
    return add(mul(sincos[0],v),mul(sincos[1],perp(v)))
def supports(v,e,sincos,pre):
    u=rotate(pre,sincos)
    pp=[]
    cr=[cross(x,u) for x in e]
    for j,(k,c) in enumerate(sector):
        if not cr[j]<0 and not cr[(j+1)%10]>0:
            pp.append(dot(v[k],u)+c)
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

def placement(v,e,a,b):
    sn=trig(box(a-300,b-300))
    a0=half-supports(v,e,sn,n0)
    a1=half-supports(v,e,sn,n1)
    trans=(a0,(2*a1-a0)/r3)
    margin=[]
    for _,u,c in rows:
        margin.append(c-supports(v,e,sn,u)-dot(u,trans))
    return margin

cells = 30
angular_cells = 100
half_bracket = F(31,25)
# The DE contact N(t0) moves by at most eps in its unit line coordinate.
# A cone with consecutive unit normals of dot product > 3/4 has half-angle
# secant below sqrt(8/7) < 107/100.
contact_error = eps
support_error = F(107,100)*contact_error
gap=3*support_error
best_root=arb(1)
best_edge=arb(1)
for j in range(27):
    for i in range(cells):
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
        predict=((1-wcen)*root_seeds[j]+wcen*root_seeds[j+1])/1000
        al,bl=predict-half_bracket,predict+half_bracket
        fs=[F_at(vertices,edges,al)-A(3*support_error),
            -F_at(vertices,edges,bl)-A(3*support_error)]
        assert min(f.lower() for f in fs)>0, (a,al,bl,fs)
        best_root=min(best_root,*(f.lower() for f in fs))
        for p in range(angular_cells):
            angle_left=al+(bl-al)*F(p,angular_cells)
            angle_right=al+(bl-al)*F(p+1,angular_cells)
            ms=placement(vertices,edges,angle_left,angle_right)
            worst=min(range(len(ms)),key=lambda k:ms[k].lower())
            assert ms[worst]>A(gap), (a,p,rows[worst][0],ms[worst],gap)
            best_edge=min(best_edge,ms[worst].lower())


r86_path = here / "eh-joint-1p5495-105p55.json"
r86 = json.loads(r86_path.read_bytes())
assert F(r86["summary"]["sigma_degrees"])==F("1.5495")
assert F(r86["summary"]["delta"])==F("0.00010555")
assert r86["summary"]["C3_arc_retained_exactly"]
assert "Ta-Z" in r86["edges"]
# Replaying those R86 duals gives the side-specific Ta-Z lower bound in
# the proof. The active K-C3 bound is irrelevant to this angle change.
assert 2*delta<F("0.0090554373129026602")
old_upper=F(r86["summary"]["candidate_area_interval_exact"][1])
area_upper=old_upper-gain_lower
assert area_upper<F("0.844093559167933607")

print(json.dumps({
    "python":sys.version, "python_flint":flint_version, "arb_precision_bits":ctx.prec,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "geometry_and_sector_source_sha256":hashlib.sha256((here/"certify_ah_Z_root_rescue.py").read_bytes()).hexdigest(),
    "R86_certificate_sha256":hashlib.sha256(r86_path.read_bytes()).hexdigest(),
    "sigma_degrees_exact":str(F("1.5495")), "AH_cut_angle_radians_exact":str(delta),
    "AH_cut_pivot_parameter":str(wp), "AH_cut_paired_test_parameter":str(wf),
    "AH_cap_height_bound_ball":str(hmax),
    "AH_cap_triangle_area_lower_exact":str(gain_lower),
    "AH_cap_triangle_area_ball":str(gain),
    "AH_graph_deficit_upper_exact":str(eps),
    "AH_graph_normal_speed_lower":str(min_eta_speed),
    "AH_graph_direct_chart_derivative_norm_upper":[str(max_xl),str(max_xt)],
    "AH_pairing_derivative_abs_upper":str(psi_max),
    "double_AH_contact_error_upper_exact":str(contact_error),
    "perturbed_compatibility_sign_slack_lower":str(best_root),
    "ideal_support_placement_slack_lower":str(best_edge),
    "support_placement_error_upper_exact":str(gap),
    "support_parameter_cells":27*cells,
    "compatibility_bracket_half_width_degrees":str(half_bracket),
    "angle_cells_per_compatibility_bracket":angular_cells,
    "directly_tested_support_rows":[label for label,_,_ in rows],
    "strict_Sprague_endpoint_Voronoi_slacks_lower":arc_signs,
    "new_area_upper_exact":str(area_upper),
    "new_area_upper_rounded_up_decimal":"0.844093559167933607",
    "passed":True,
},indent=2))
