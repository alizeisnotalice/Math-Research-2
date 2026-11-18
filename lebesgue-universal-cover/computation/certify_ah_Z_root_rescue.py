"""Exact support cones and a new compatible placement for limiting double AH hits."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

source=Path(__file__).with_name("certify_gibbs_area_refined.py")
raw=source.read_bytes()
parts=raw.decode().split("width = F(1,cells)")
assert len(parts)==2
sys.argv=[str(source),"1.5495","2","1"]
d={"__file__":str(source)}
exec(compile(parts[0],str(source),"exec"),d)
I,add,sub,mul,dot,cross,perp=(d[k] for k in ("I","add","sub","mul","dot","cross","perp"))
trig,pi=d["trig"],d["pi_machin"]
half,r3=d["half"],d["r3"]
n=d["n"]
l,z,mid=d["l0"],d["z"],d["mid"]

def squares(p):
    return dot(p,p)

def unit_degrees(q):
    q=I(q)
    return tuple(I(*trig.trig(trig.mul(pi,(q.lo/180,q.hi/180)),k)) for k in ("cos","sin"))

def rotation(p,co,si):
    return add(mul(co,p),mul(si,perp(p)))

def normalize(p):
    return mul(1/squares(p).sqrt(),p)

def qa(v):
    return (half*(r3*v[1]-v[0]),half*(r3*v[0]+v[1]))

v=[l,qa(l),qa(z),mid,z]
edges={}
for i in range(5):
    for j in [(i+1)%5,(i-1)%5]:
        edges[i,j]=sub(v[i],v[j])
        check=squares(edges[i,j])
        assert check.lo<=1<=check.hi
for i in range(5):
    for j in range(i+1,5):
        if (i,j) not in edges:
            assert squares(sub(v[i],v[j])).hi<1

# Every adjacent pair in this ordered normal list must be in strict
# CCW order, with the first five in the upper half-plane.
order=[(2,1),(2,3),(4,3),(4,0),(1,0),(1,2),(3,2),(3,4),(0,4),(0,1)]
for idx,edge in enumerate(order):
    vv=edges[edge]
    assert (vv[1].lo>0 if idx<5 else vv[1].hi<0)
    ww=edges[order[(idx+1)%10]]
    assert cross(vv,ww).lo>0 and dot(vv,ww).lo>0
sector_formula=[(2,0),(3,1),(4,0),(0,1),(1,0),(2,1),(3,0),(4,1),(0,0),(1,1)]

def support_at(direction):
    possible=[]
    for i,(center,on_arc) in enumerate(sector_formula):
        aa=edges[order[i]]
        bb=edges[order[(i+1)%10]]
        if cross(aa,direction).hi>=0 and cross(direction,bb).hi>=0:
            # Edges are tangent to two unit disks at the marked vertex.
            # An opposite direction has the exact complementary support value.
            for edge_label in (order[i],order[(i+1)%10]):
                assert edge_label[1 if on_arc else 0]==center
            possible.append(dot(v[center],direction)+on_arc)
    assert possible
    return I(min(x.lo for x in possible),max(x.hi for x in possible))

def F_compat(phi):
    terms=[support_at(unit_degrees(I(phi)+60*j)) for j in range(3)]
    return terms[0]-terms[1]+terms[2]-half

lo_phi,hi_phi=F("29.997"),F("30.003")
left,right=F_compat(lo_phi),F_compat(hi_phi)
assert left.hi<-F(23,10**7)<0<F(23,10**7)<right.lo
pre=I(lo_phi-240,hi_phi-240)
rot_cos,rot_sin=unit_degrees(pre)
premap=lambda x:rotation(x,rot_cos,rot_sin)
a0=half-support_at(premap(n[0]))
a1=half-support_at(premap(n[1]))
translation=a0,(2*a1-a0)/r3

# Build the inscribed polygon using the exact R86 EH depth and 4
# radial interpolants per short arc. Its six hexagon faces are exact.
u=d["u"]
eh_delta=F(10555,10**8)
co,si=tuple(I(*trig.trig((eh_delta,eh_delta),k)) for k in ("cos","sin"))
c3=d["c3"]
u1=add(c3,rotation(sub(u,c3),co,-si))
def p_curve(w):
    return d["pair"](add(z,mul(I(w),d["direction"])))
def arcs(a,b,c):
    aa,bb=sub(a,c),sub(b,c)
    assert cross(aa,bb).lo>0 and dot(aa,bb).lo>0
    return [add(c,normalize(add(mul(I(1-F(i,4)),aa),mul(I(F(i,4)),bb)))) for i in range(1,4)]+[b]
tp=p_curve(F(1))
q,y,t,e3,c,g,p,f3,kk,h,te,ip,b3=(d[k] for k in ("q","y","t","e3","c","g","p","f3","k","h","t_eh","ip","b3"))
def edge_list(label,point,tail):
    rows=[]
    for j,after in enumerate(tail):
        normal=normalize(mul(-1,perp(sub(after,point))))
        rows.append((f"{label}-{j}",normal,dot(normal,point)))
        point=after
    return rows

rows=(edge_list("Q-Y",q,arcs(q,y,e3))
      +edge_list("Y-Ta",y,arcs(y,t,c))
      +edge_list("Ta-Z",t,[z])
      +edge_list("Z-Tprime",z,[p_curve(F(j,4)) for j in range(1,5)])
      +edge_list("Tprime-P",tp,arcs(tp,p,g))
      +edge_list("H-U1",h,arcs(h,u1,c3))
      +edge_list("U1-Te",u1,[te])
      +edge_list("Te-I",te,arcs(te,ip,b3))
      +edge_list("I-E3",ip,[e3])
      +edge_list("G-K",g,arcs(g,kk,f3))
      +edge_list("K-C3",kk,[c3]))
assert len(rows)==32
margins={}
for label,normal,intercept in rows:
    support=support_at(premap(normal))+dot(normal,translation)
    margin=intercept-support
    assert margin.lo>F(1,10000),(label,margin.out())
    margins[label]=str(margin.lo)

print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "coordinate_source_sha256":hashlib.sha256(raw).hexdigest(),
    "sigma_degrees_exact":str(d["sigma"]),
    "baseline_EH_depth_radians_exact":str(eh_delta),
    "compatibility_range_exact_degrees":[str(lo_phi),str(hi_phi)],
    "F_left_interval":left.out(),"F_right_interval":right.out(),
    "vertex_or_antipodal_arc_sectors":[[list(order[i]),list(order[(i+1)%10]),*formula]
                                         for i,formula in enumerate(sector_formula)],
    "new_theta_expression":"240 degrees - phi",
    "nonhex_edges_strict_margins":margins,
    "minimum_nonhex_edge_margin_lower":str(min(map(F,margins.values()))),
    "passed":True,
},indent=2))
